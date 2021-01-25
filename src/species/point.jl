#==============================================================================#
# Point-like particle
#==============================================================================#

"""
	mutable struct Point <: AbstractSpecies

Species corresponding to a single point-like particle composed of a single monomer type.
"""
mutable struct Point <: AbstractSpecies
	mids     :: Vector{Int}
	Nref     :: Float64

	# SCF fields
	Q        :: Float64
	density  :: Dict{Int,FieldGrid{Float64}}
	system   :: Option{FieldSystem}

	function Point(mon::Monomer)
		return new([mon.id], mon.vol, 0.0, Dict(), nothing)
	end
end

#==============================================================================#

Base.show(io::IO, point::Point) = @printf(io, "Point(mid = %d)", point.mids[1])

monomer_fraction(point::Point, mid::Integer) = hasmonomer(point, mid) ? 1.0 : 0.0

function setup!(point::Point, sys::FieldSystem)
	point.system = sys
	point.Q = 0.0
	point.density[point.mids[1]] = zeros(Float64, sys.dims)
end

function density!(point::Point)
	@assert !isnothing(point.system)
	sys = point.system
	mon = sys.monomers[point.mids[1]]
	omega = sys.fields[mon.id]

	# Partition function via Boltzmann weight
	Q = 0.0
	@simd for i in eachindex(omega)
		@inbounds Q += exp(-mon.vol * omega[i])
	end
	point.Q = Q / ngrid(sys) # Normalize by number of grid points

	# Update point density field
	@. point.density[mon.id] = exp(-mon.vol * omega) / point.Q

	return nothing
end

function scfstress(point::Point)
	@assert !isnothing(point.system)
	return zeros(nparams(sys.cell))
end