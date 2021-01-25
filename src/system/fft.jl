"""
	struct FFTHolder{FP,IP}

A container for pre-allocated FFT grids and plans for a provided grid-size.
"""
struct FFTHolder{FP,IP}
    dims   :: NTuple{3,Int}
    rgrids :: Vector{FieldGrid{Float64}}
    kgrids :: Vector{FieldGrid{Complex{Float64}}}

    FT     :: FP
    iFT    :: IP
	# FT     :: FFTW.rFFTWPlan{Float64,-1,false,3}
    # iFT    :: AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.rFFTWPlan{Complex{Float64},1,false,3},Float64}
end

function FFTHolder(dims::NTuple{3,<:Integer}; nthreads::Integer = -1)
	gr = zeros(Float64, dims)
	gk = zeros(Complex{Float64}, ksize(dims))

	TGR = typeof(gr)
	TGK = typeof(gk)

	# Setup multi-threads and FFT plan
	if nthreads > 0
		FFTW.set_num_threads(nthreads)
	else
		nthreads = determine_nthreads(prod(dims))
		FFTW.set_num_threads(nthreads)
	end
	FT = plan_rfft(gr; flags = FFTW.MEASURE)
	iFT = plan_irfft(gk, dims[1]; flags = FFTW.MEASURE) #inv(FT)
	plan = FFTHolder(dims, TGR[], TGK[], FT, iFT)

	return plan
end

FFTHolder(Nx::Integer, Ny::Integer, Nz::Integer; kwargs...) = FFTHolder((Nx, Ny, Nz); kwargs...)

#==============================================================================#

Base.show(io::IO, plan::FFTHolder) = @printf(io, "%s(dims = %s)", typeof(plan), plan.dims)

Base.empty!(plan::FFTHolder) = (empty!(plan.rgrids); empty!(plan.kgrids);)

nrgrids(plan::FFTHolder) = length(plan.rgrids)
nkgrids(plan::FFTHolder) = length(plan.kgrids)

add_rgrid!(plan::FFTHolder) = push!(plan.rgrids, zeros(Float64, plan.dims))
add_kgrid!(plan::FFTHolder) = push!(plan.kgrids, zeros(Complex{Float64}, ksize(plan.dims)))

function rgrid(plan::FFTHolder, idx::Integer)
	while nrgrids(plan) < idx
		add_rgrid!(plan)
	end
	return plan.rgrids[idx]
end

function kgrid(plan::FFTHolder, idx::Integer)
	while nkgrids(plan) < idx
		add_kgrid!(plan)
	end
	return plan.kgrids[idx]
end