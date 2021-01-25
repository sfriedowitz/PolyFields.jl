#==============================================================================#
# System
#==============================================================================#

"""
	@enum Ensemble::UInt8 Canonical Grand

Statistical ensemble for a field-theoretic system.
"""
@enum Ensemble::UInt8 begin
	Canonical
    Grand
end

"""
	mutable struct FieldSystem <: AbstractSystem

A system containing data for a field-theoretic simulation,
represented as a set of species, field grids, and a corresponding unit cell.
Contains all species, field, density, and interaction information.
"""
mutable struct FieldSystem <: AbstractSystem
	# Grid and cell
	dims          :: NTuple{3,Int}
	ensemble      :: Ensemble 
	cell          :: UnitCell

	# MDE and FFT
	fftplan       :: FFTHolder
	solver        :: PseudoSpectralSolver

	# Field grids for each monomer type
	monomers      :: Dict{Int,Monomer}
	monomer_fracs :: Dict{Int,Float64}
	fields        :: Dict{Int,FieldGrid{Float64}}
	potentials    :: Dict{Int,FieldGrid{Float64}}
	residuals     :: Dict{Int,FieldGrid{Float64}}
	density       :: Dict{Int,FieldGrid{Float64}}

	# Interactions and species
	#constraints  :: Vector{AbstractConstraint}
	interactions  :: Vector{AbstractInteraction}
	species       :: Vector{AbstractSpecies}
	phi_species   :: Vector{Float64}
	mu_species    :: Vector{Float64}

	# Constraints
	compress      :: Compressibility
end

function FieldSystem(dims::NTuple{3,<:Integer}, cell::UnitCell;
	monomers = Monomer[], ensemble::Ensemble = Canonical, 
	compress = Compressibility(0.0), mde::Symbol = :RK2, nthreads::Integer = -1)
	# Create FFT and MDE solver
	fftplan = FFTHolder(dims; nthreads = nthreads)
	solver = PseudoSpectralSolver(dims; method = mde)

	sys = FieldSystem(dims, ensemble, cell, fftplan, solver,
		Dict(), Dict(), Dict(), Dict(), Dict(), Dict(), [], [], [], [], 
		compress
	)

	# Allocate to system size
	setup!(cell, sys)
	setup!(compress, sys)

	# Add initial monomers
	for mon in monomers
		add_monomer!(sys, mon)
	end

	return sys
end

FieldSystem(nx::Integer, ny::Integer, nz::Integer, cell::UnitCell; kwargs...) = FieldSystem((nx, ny, nz), cell; kwargs...)

#==============================================================================#

function Base.show(io::IO, sys::FieldSystem)
	@printf(io, "FieldSystem(dims = %s, ensemble = %s, %d monomers, %d species)", 
		sys.dims, sys.ensemble, nmonomers(sys), nspecies(sys))
end

volume(sys::FieldSystem) = sys.cell.volume

ngrid(sys::FieldSystem) = prod(sys.dims)

nspecies(sys::FieldSystem) = length(sys.species)

nmonomers(sys::FieldSystem) = length(sys.monomers)

hasmonomer(sys::FieldSystem, mid::Integer) = haskey(sys.monomers, mid)

#==============================================================================#

"""
	add_cell!(sys, cell)

Provide a new cell for the system.
"""
function add_cell!(sys::FieldSystem, cell::UnitCell)
	sys.cell = cell
	setup!(cell, sys)
	return nothing
end

"""
	add_monomer!(sys, mon)

Add a `Monomer` to the `FieldSystem`.
"""
function add_monomer!(sys::FieldSystem, mon::Monomer)
	if hasmonomer(sys, mon.id)
		@warn "Replacing monomer with id = $(mon.id) in FieldSystem."
	end

	# Add to system, add new fields
	sys.monomers[mon.id] = mon
	sys.monomer_fracs[mon.id] = 0.0

	sys.fields[mon.id] = zeros(Float64, sys.dims)
	sys.potentials[mon.id] = zeros(Float64, sys.dims)
	sys.residuals[mon.id] = zeros(Float64, sys.dims)
	sys.density[mon.id] = zeros(Float64, sys.dims)

	return nothing
end

"""
	add_species!(sys, species, phi/mu)

Add an `AbstractSpecies` to the `FieldSystem`.
For `Canonical` ensemble systems, the third argument is treated as the species volume fraction.
For `Grand` ensemble systems, this argument is the species chemical potential.
"""
function add_species!(sys::FieldSystem, species::AbstractSpecies, val::Real)
	# Check to make sure species has valid monomers
	for mid in species.mids
		if !hasmonomer(sys, mid)
			error("Species contains monomer with id = $(mid) not present in system.")
		end
	end

	setup!(species, sys)
	push!(sys.species, species)

	# Add value based on ensemble of simulation
	if sys.ensemble == Canonical
		if val <= 0.0
			error("Bulk volume fraction must be 0 < phi <= 1.")
		end
		push!(sys.phi_species, val)
		push!(sys.mu_species, 0.0)
	elseif sys.ensemble == Grand
		push!(sys.phi_species, 0.0)
		push!(sys.mu_species, val)
	end

	return nothing
end

"""
	add_interaction!(sys, itx)

Add an `AbstractInteraction` to the `FieldSystem`.
"""
function add_interaction!(sys::FieldSystem, itx::AbstractInteraction)
	# Check to make sure has valid monomers
	for mid in itx.mids
        if !hasmonomer(sys, mid)
            error("Interaction contains monomer with id = $(mid) not present in system.")
        end
    end

	setup!(itx, sys)
	push!(sys.interactions, itx)
	return nothing
end

"""
	validate(sys)

Determine if the system is properly initialized before a calculation.
Checks if the following requirements are met:

* All fields are initialized with non-zero mean or standard deviation
* All monomer types are found in at least one species
* Bulk fractions of all monomers sums to unity (canonical ensemble)
* Bulk charge fractions of all monomers sums to zero
"""
function validate(sys::FieldSystem)
	# Check if fields are initialized
	fields_init = false
	for omega in values(sys.fields)
		if mean(omega) != 0.0 || std(omega) != 0.0
			fields_init = true
			break
		end
	end
	if !fields_init
		error("Fields must be initialized before SCFT simulation.")
	end

	# Check bulk volume fraction constraint
	if sys.ensemble == Canonical
		if !(sum(sys.phi_species) ≈ 1.0)
			error("Species volume fractions do not sum to unity.")
		end
	end

	# Check bulk charge neutrality
	bulk_charge = 0.0
	monomer_fractions!(sys)
	for (mid, mon) in sys.monomers
		bulk_charge += sys.monomer_fracs[mid] * mon.charge / mon.vol		
	end
	if !(bulk_charge ≈ 0.0)
		error("Bulk system not charge neutral.")
	end

	# Check if any monomer types do not have associated species
	for mid in keys(sys.monomers)
		has_species = false
		for species in sys.species
			if hasmonomer(species, mid)
				has_species = true
				break
			end
		end
		if !has_species
			error("Monomer with id = $(mid) not found in any species.")
		end
	end

	return nothing
end

#==============================================================================#
# Methods for system density, potentials, and constraints
#==============================================================================#

"""
	muphi!(sys)

Update species volume fractions and chemical potentials.
"""
function muphi!(sys::FieldSystem)
	if sys.ensemble == Canonical
		for (isp, species) in enumerate(sys.species)
			sys.mu_species[isp] = log(sys.phi_species[isp] / species.Q)
		end
	elseif sys.ensemble == Grand
		for (isp, species) in enumerate(sys.species)
			sys.phi_species[isp] = species.Q * exp(sys.mu_species[isp])
		end
	end
	return nothing
end

"""
	density!(sys)

Compute density fields for all monomers in the system.
"""
function density!(sys::FieldSystem)
	# Zero density fields
	for rho in values(sys.density)
		rho .= 0.0
	end

	# Loop over all species, update density on block/type level
	for species in sys.species
		density!(species)
	end

	# Update all bulk fractions/chemical potentials
	muphi!(sys)

	# Pool the updated density by type at the system level
	for (isp, species) in enumerate(sys.species)
		phi = sys.phi_species[isp]
		for mid in species.mids
			@. sys.density[mid] += phi * species.density[mid]
		end
	end

	return nothing
end

"""
	potentials!(sys)

Compute differentials of the excess interactions with respect to each density field.
"""
function potentials!(sys::FieldSystem)
	# Zero current potentials
	for pot in values(sys.potentials)
		pot .= 0.0
	end

	# Add excess interactions
	for (mid, pot) in sys.potentials
		for itx in sys.interactions
			potential!(itx, mid, pot)
		end
	end

	return nothing
end

"""
    residuals!(sys)

Compute the gradient of the system Hamiltonian with respect to each density field.
"""
function residuals!(sys::FieldSystem)
	# Run updates on all field grids
	density!(sys)
	potentials!(sys)

	# Update compressibility field
	update!(sys.compress)
	eta = sys.compress.field

    # Make residuals from fields and calculated potentials/eta field
	for (mid, res) in sys.residuals
		@. res = sys.potentials[mid] + eta - sys.fields[mid]
	end

	# Mean-subtract the system fields and residuals
	if sys.ensemble == Canonical
		for mid in keys(sys.fields)
			sys.fields[mid] .-= mean(sys.fields[mid])
			sys.residuals[mid] .-= mean(sys.residuals[mid])
		end
	end
    
	return nothing
end

"""
	uniform_fields!(sys)

Set the value of all chemical potential fields equal to the summation of interaction potentials, 
corresponding to the mean-field limit with vanishing constraint fields.
"""
function uniform_fields!(sys::FieldSystem)
	potentials!(sys)
	for (mid, omega) in sys.fields
		@. omega = sys.potentials[mid]
		if sys.ensemble == Canonical; 
			omega .-= mean(omega); 
		end
	end
	return nothing	
end

"""
	monomer_fractions!(sys)

Update the bulk fractions for each monomer present in the system.
"""
function monomer_fractions!(sys::FieldSystem)
	bulk = sys.monomer_fracs
	for mid in keys(bulk); bulk[mid] = 0.0; end

	for (isp, species) in enumerate(sys.species)
		phi = sys.phi_species[isp]
		for mid in species.mids
			bulk[mid] += phi * monomer_fraction(species, mid)
		end
	end
	return bulk
end

#==============================================================================#

"""
	scfstress(sys)

Compute the free energy variation of the system with respect to each cell parameter.
"""
function scfstress(sys::FieldSystem)
	stress = zeros(nparams(sys.cell))
	muphi!(sys)

	# Add stress from each contributing species
	for (isp, species) in enumerate(sys.species)
		phi = sys.phi_species[isp]
		stress .-= phi * scfstress(species)
	end

	return stress
end

"""
	free_energy(sys)

Return the intensive Helmholtz free energy for the interacting, inhomogeneous system.
"""
function free_energy(sys::FieldSystem)
	# Add translational entropy terms
	ftrans = 0.0
	for (isp, species) in enumerate(sys.species)
		phi = sys.phi_species[isp]
		mu = sys.phi_species[isp]
		if phi > 1e-8
			ftrans += phi * (mu - 1.0) / species.Nref
		end
	end

	# Add field/density pair
	ffield = 0.0
	for (mid, omega) in sys.fields
		rho = sys.density[mid]
		@simd for i in eachindex(rho)
			@inbounds ffield -= rho[i] * omega[i]
		end
	end
	ffield /= ngrid(sys)

	# Add excess interaction terms
	fexc = 0.0
	for itx in sys.interactions
		fexc += energy(itx)
	end

	# All pieces together
	fhelm = ftrans + ffield + fexc

	return fhelm
end

"""
	free_energy_bulk(sys)

Return the intensive Helmholtz free energy for a system at the specified bulk composition.
"""
function free_energy_bulk(sys::FieldSystem)
	fhelm = 0.0

	# Add translational entropy terms
	for (isp, species) in enumerate(sys.species)
		phi = sys.phi_species[isp]
		if phi > 1e-8
			fhelm += phi * (log(phi) - 1.0) / species.Nref
		end
	end

	# Add the excess interaction energy
	monomer_fractions!(sys)
	for itx in sys.interactions
		fhelm += energy_bulk(itx)
	end

	return fhelm
end

"""
    field_error(sys)

Compute the weighted error based on the current field residuals.

Follows the definition in the Matsen (2009) Anderson mixing papers:
``
{\\rm err} = [∑_{α,i} {\\rm res}_{α,i}^2 / ∑_{α,i} ω_{α,i}^2]^{1/2}
``
"""
function field_error(sys::FieldSystem)
    # Residual weights
    rsum = 0.0
    for (mid, res) in sys.residuals
        @simd for i in eachindex(res)
            @inbounds rsum += res[i]^2
        end
    end
    # Field weights
    wsum = eps(Float64)
    for (mid, omega) in sys.fields
        @simd for i in eachindex(omega)
            @inbounds wsum += omega[i]^2
        end
    end
    return sqrt(rsum) / sqrt(wsum)
end