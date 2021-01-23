#==============================================================================#
# Point-like particle
#==============================================================================#

"""
	mutable struct Point <: AbstractSpecies

Species corresponding to a single point-like particle composed of a single monomer type.
"""
mutable struct Point <: AbstractSpecies
	mids    :: Vector{Int}
	Nref    :: Float64

	# SCF fields
	Q       :: Float64
	density :: Dict{Int,FieldGrid{Float64}}
	system  :: Option{FieldSystem}

	function Point(mon::Monomer)
		return new(name, [mon.id], mon.size, 0.0, Dict(), nothing)
	end
end

#==============================================================================#
# Methods
#==============================================================================#

Base.show(io::IO, species::Point) = @printf(io, "Point(mid = %d)", species.mids[1])

monomer_fraction(point::Point, mid::Integer) = has_monomer(point, mid) ? 1.0 : 0.0

function setup!(point::Point, sys::FieldSystem)
	point.system = sys
	point.Q = 0.0
	point.density[point.mids[1]] = zeros(Float64, sys.npw)
	return nothing
end

function density!(species::Point)
	@assert !isnothing(point.system)
	sys = point.system
	mon = sys.monomers[point.mids[1]]
	omega = sys.fields[mon.id]

	# Partition function via Boltzmann weight
	Q = 0.0
	@simd for i in eachindex(omega)
		@inbounds Q += exp(-mon.size * omega[i])
	end
	point.Q = Q / num_grid(sys) # Normalize by number of grid points

	# Update species density field
	@. point.density[mon.id] = exp(-mon.size * omega) / point.Q

	return nothing
end

function scf_stress(point::Point)
	return zeros(num_cell_params(point.system.cell))
end