"""
	mutable struct Multiblock <: AbstractSpecies

A species representing a linear polymer chain with variable chain length and number of blocks.
"""
mutable struct Multiblock <: AbstractSpecies
	monomers      :: Vector{Monomer}
	mids          :: Vector{Int}
	mid_block     :: Vector{Int}	       # Block idx  --> monomer id
	mid_map       :: Dict{Int,Vector{Int}} # Monomer id --> block idx, blocks which contain a given mid

	# Chain structure
	Ns            :: Int
    N             :: Int
    N_block       :: Vector{Int}
    Nref          :: Float64
    Nref_block    :: Vector{Float64}
	f_block       :: Vector{Float64}
	b_block       :: Vector{Float64}
	ds_block      :: Vector{Float64}
	begin_block   :: Vector{Int}

	# SCF fields
	Q             :: Float64
	q             :: Vector{FieldGrid{Float64}}
	qc            :: Vector{FieldGrid{Float64}}
	density       :: Dict{Int,FieldGrid{Float64}} # Stores density per monomer type
	density_block :: Dict{Int,FieldGrid{Float64}} # Stores density per block idx
	system        :: Option{FieldSystem}
end

function Multiblock(mon_block::AbstractVector{Monomer}, N_block::AbstractVector{<:Integer}, 
    b_block::AbstractVector{<:Real}, Ns::Integer = 100)
    @assert length(mon_block) == length(N_block) == length(b_block)

    # Setup monomer info and block-type mappings
    monomers = sort(unique(mon_block))
	mids = [mon.id for mon in monomers]

    mid_block = [mon.id for mon in mon_block]
    mid_map = Dict()
    for (iblk, mid) in enumerate(mid_block)
        if !haskey(mid_map, mid)
            mid_map[mid] = []
        end
        push!(mid_map[mid], iblk)
    end
    
    # Determine contour locations for each block
    N = sum(N_block)
    Nref_block = [mon_block[iblk].vol * N_block[iblk] for iblk = 1:length(mon_block)]
    Nref, begin_block, ds_block = chain_grid(Nref_block, Ns) # Use the Nref block lengths to discretize solver
    f_block = [Nr/Nref for Nr in Nref_block]

    # Check if values were adjusted
    if Ns != begin_block[end]-1
        @warn "Non-integer chain discretization. Adjusting `ds` and `Ns` appropriately."
        Ns = begin_block[end]-1
    end

    return Multiblock(monomers, mids, mid_block, mid_map,
        Ns, N, N_block, Nref, Nref_block, f_block, b_block, ds_block, begin_block,
        0.0, [], [], Dict(), Dict(), nothing
    )
end

function Homopolymer(mon::Monomer, N::Integer, b::Real, Ns::Integer = 100)
    return Multiblock([mon], [N], [b], Ns)
end

function Diblock(amon::Monomer, bmon::Monomer, N::Integer, f::Real, b1::Real, b2::Real, Ns::Integer = 100)
    N_block = trunc.(Int, [f*N, (1-f)*N])
    if sum(N_block) != N
        error("Non-integer block lengths for N = $N, f = $f.")
    end
    return Multiblock([amon, bmon], N_block, [b1, b2], Ns)
end

#==============================================================================#

function Base.show(io::IO, chain::Multiblock)
	if nblocks(chain) == 1
        @printf(io, "Homopolymer(mid = %d, N = %d)", chain.mids[1], chain.N)
    elseif nblocks(chain) == 2
        @printf(io, "Diblock(mids = %s, N = %d, f = [%.3f, %.3f])", chain.mids, chain.N, chain.f_block[1], chain.f_block[2])
    else
        @printf(io, "Multiblock(%d blocks, %d monomers, N = %d)", num_monomers(chain), nblocks(chain), chain.N)
    end
end

nblocks(chain::Multiblock) = length(chain.N_block)

function monomer_fraction(chain::Multiblock, mid::Integer)
	# Check if mid present
	if !hasmonomer(chain, mid); return 0.0; end
	# It has the mid, so calculate fraction
	f = 0.0
	for iblk in chain.mid_map[mid]
		f += chain.f_block[iblk]
	end
	return f
end

function setup!(chain::Multiblock, sys::FieldSystem)
    chain.system = sys

	# Setup propagators
	chain.Q = 0.0
	chain.q = [zeros(Float64, sys.dims) for n = 1:chain.Ns+1]
	chain.qc = [zeros(Float64, sys.dims) for n = 1:chain.Ns+1]

	# Setup density grids
	for mid in chain.mids
		chain.density[mid] = zeros(Float64, sys.dims)
	end
	for iblk = 1:nblocks(chain)
		chain.density_block[iblk] = zeros(Float64, sys.dims)
	end

	return nothing
end

#==============================================================================#

function density!(chain::Multiblock)
	@assert !isnothing(chain.system)
	sys = chain.system
    plan = sys.fftplan

	# Necessary class fields
	q, qc = chain.q, chain.qc
	solver = sys.solver
	ksq = sys.cell.ksq

	# Solve q by integrating forward along chain
	q[1] .= 1.0
	for iblk = 1:nblocks(chain)
		b = chain.b_block[iblk]
		ds = chain.ds_block[iblk]

		mid = chain.mid_block[iblk]
		mon = sys.monomers[mid]
		omega = sys.fields[mid]

		update!(solver, omega, ksq, mon.vol, b, ds)
		for s = chain.begin_block[iblk] : chain.begin_block[iblk+1]-1
			propagate!(solver, plan, q[s], q[s+1])
		end
	end

	# Solve qc by integrating reverse along the chain
    # If a homopolymer, short-cut the integration and copy the propagator in reverse
    if nblocks(chain) == 1
        for (s, qi) in enumerate(Iterators.reverse(q))
            @inbounds qc[s] .= qi
        end
    else
    	qc[end] .= 1.0
    	for iblk = nblocks(chain):-1:1
    		b = chain.b_block[iblk]
    		ds = chain.ds_block[iblk]

			mid = chain.mid_block[iblk]
			mon = sys.monomers[mid]
			omega = sys.fields[mid]

    		update!(solver, omega, ksq, mon.vol, b, ds)
    		for s = chain.begin_block[iblk+1] : -1 : chain.begin_block[iblk]+1
    			propagate!(solver, plan, qc[s], qc[s-1])
    		end
    	end
    end

	# Partition function
	chain.Q = sum(chain.q[end]) / ngrid(sys)

	# Integrate density w/ Simpson's rule
	# Reset density grids
	for rho in values(chain.density); rho .= 0.0; end
	for rho in values(chain.density_block); rho .= 0.0; end

	for iblk = 1:nblocks(chain)
		# Block start and end indices
		blk_bgn = chain.begin_block[iblk]
		blk_end = chain.begin_block[iblk+1]

		ds = chain.ds_block[iblk]
		rho = chain.density_block[iblk]

		# First and last
		@simd for i in eachindex(rho)
			@inbounds rho[i] += q[blk_bgn][i] * qc[blk_bgn][i]
		end
		@simd for i in eachindex(rho)
			@inbounds rho[i] += q[blk_end][i] * qc[blk_end][i]
		end

		# Odd indices
		for s = blk_bgn+1:2:blk_end-1
            qs, qcs = q[s], qc[s]
			@simd for i in eachindex(rho)
				@inbounds rho[i] += 4.0 * qs[i] * qcs[i]
			end
		end

		# Even indices
		for s = blk_bgn+2:2:blk_end-2
            qs, qcs = q[s], qc[s]
			@simd for i in eachindex(rho)
				@inbounds rho[i] += 2.0 * qs[i] * qcs[i]
			end
		end

		# Normalize by ds/3 for Simpson's rule
		rho .*= (ds / 3.0)
	end

    # Normalize by 1/(Q*N)
    for rho in values(chain.density_block)
        rho ./= chain.Q * chain.Nref
    end

	# Add each block density to appropriate monomer type
	for mid in chain.mids
		for iblk in chain.mid_map[mid]
			chain.density[mid] .+= chain.density_block[iblk]
		end
	end

	return nothing
end

function scfstress(chain::Multiblock)
	@assert !isnothing(chain.system)
	sys = chain.system
	cell = sys.cell
	plan = sys.fftplan
	q, qc = chain.q, chain.qc

    # Get temp storage grids from the FFT helper
    qk1, qk2, qtmp = kgrid(plan,1), kgrid(plan,2), kgrid(plan,3)

	# Calculate variation in partition function by looping over all chain contour segments
	dQ = zeros(nparams(cell))
	for iblk = 1:nblocks(chain)
		# Block start and end indices
		blk_bgn = chain.begin_block[iblk]
		blk_end = chain.begin_block[iblk+1]

		b = chain.b_block[iblk]
		ds0 = chain.ds_block[iblk]

		# Block start and end indices
		for s = blk_bgn:blk_end
			# Adjust prefactor ds for Simpson's rule quadrature
			ds = ds0
			if s != blk_bgn && s!= blk_end
				if s % 2 == 0
					ds *= 4.0
				else
					ds *= 2.0
				end
			end

            # FFT the propagators to k-space grids
            qs = q[s]; qcs = qc[s]
            mul!(qk1, plan.FT, qs)
            mul!(qk2, plan.FT, qcs)

            # Apply dksq for each parameter
            for k = 1:nparams(cell)
                dksq = cell.dksq[k]
                @. qtmp = qk1 * dksq

                # Perform summation over kgrid propagators
                # Since we use half-space real FFT, the complex conjugate components do not cancel
                # We obtain the full result by doubling all the real components
                #   of the k != 0 waves
                qint = real(qk2[1]) * real(qtmp[1])
                @simd for i = 2:length(qk2)
                    @inbounds qint += 2 * real(qk2[i]) * real(qtmp[i])
                end

                dQ[k] -= b^2 * ds * qint / 6.0 / ngrid(sys)^2
            end

		end
	end

	# Normalize stress after completion
	dQ ./= 3.0 # Simpson's rule factor
	dQ ./= chain.Q * chain.Nref # Stress formula normalization

    return dQ
end