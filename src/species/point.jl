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
		return new([mon.id], mon.vol, 0.0, Dict())
	end
end

#==============================================================================#
# Methods
#==============================================================================#

Base.show(io::IO, species::Point) = @printf(io, "Point(mid = %d)", species.mids[1])

monomer_fraction(species::Point, mid::Integer) = hasmonomer(species, mid) ? 1.0 : 0.0

function setup!(species::Point, sys::FieldSystem)
	species.system = sys
	species.Q = 0.0
	species.density[species.mid[1]] = zeros(Float64, dims)
end

function density!(species::Point)
	@assert !isnothing(species.system)
	sys = species.system
	mon = sys.monomers[species.mids[1]]

	# Partition function via Boltzmann weight
	Q = 0.0
	@simd for i in eachindex(omega)
		@inbounds Q += exp(-mon.vol * omega[i])
	end
	species.Q = Q / ngrid(sys) # Normalize by number of grid points

	# Update species density field
	@. species.density[mon.id] = exp(-mon.vol * omega) / point.Q

	return nothing
end

scfstress(species::Point, cell::Cell) = zeros(nparams(cell))