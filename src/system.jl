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
	FieldSystem(npw, cell; ensemble, monomers, mde, nthreads)

A system containing data for a field-theoretic calculation.
Contains all species, field, density, and monomer type information.
"""
mutable struct FieldSystem
	# Grid and cell
	npw          :: NTuple{3,Int}
	ensemble     :: Ensemble 

	# Cell and stress
	cell         :: Cell
	domain       :: Bool
	stress       :: Vector{Float64}

	# FFT and MDE helpers
	fft          :: FFTBuddy
	solver       :: PseudoSpectralSolver

	# Field grids for each monomer type
	monomers     :: Dict{Int,Monomer}
	fields       :: Dict{Int,FieldGrid{Float64}}
	potentials   :: Dict{Int,FieldGrid{Float64}}
	residuals    :: Dict{Int,FieldGrid{Float64}}
	density      :: Dict{Int,FieldGrid{Float64}}
	density_bulk :: Dict{Int,Float64}

	# Constraints
	constraints  :: Vector{AbstractConstraint}
	elec         :: Electroneutrality
	eta          :: FieldGrid{Float64}
	kappa        :: FieldGrid{Float64}
	density_sum  :: FieldGrid{Float64}
	charge_sum   :: FieldGrid{Float64}

	# Species, interactions, constraints
	species      :: Vector{AbstractSpecies}
	phi_species  :: Vector{Float64}
	mu_species   :: Vector{Float64}
	interactions :: Vector{AbstractInteraction}
end

#==============================================================================#
# Constructors
#==============================================================================#

function FieldSystem(npw::NTuple{3,<:Integer}, cell::Cell; 
	domain::Bool = true, 
	ensemble::Ensemble = Canonical, 
	monomers::AbstractVector{Monomer} = Monomer[],
	comp::Compressibility = Compressibility(),
	mde::AbstractString = "RK2", 
	nthreads::Integer = -1
)

	# Setup k-vector grids in the cell
	setup_kgrids!(cell, npw)
	stress = zeros(num_cell_params(cell))

	# Create FFT and MDE solver
	fft = FFTBuddy(npw; ngrids = 0, nthreads = nthreads)
	solver = PseudoSpectralSolver(npw; method = mde)

	# Create system
	sys = FieldSystem(npw, ensemble, cell, domain, stress, fft, solver,
		Dict(), Dict(), Dict(), Dict(), Dict(), Dict(),
		comp, zeros(npw), zeros(npw),
		[], [], [], []
	)

	# Add provided monomers
	for mon in monomers; add_monomer!(sys, mon); end

	return sys
end

FieldSystem(Nx::Integer, Ny::Integer, Nz::Integer, cell::Cell; kwargs...) = 
	FieldSystem((Nx, Ny, Nz), cell; kwargs...)

#==============================================================================#
# Methods
#==============================================================================#

function Base.show(io::IO, sys::FieldSystem)
	@printf(io, "FieldSystem(npw = %s, ensemble = %s, %d monomers, %d species)", 
		sys.npw, sys.ensemble, num_monomers(sys), num_species(sys))
end

num_grid(sys::FieldSystem) = prod(sys.npw)
num_monomers(sys::FieldSystem) = length(sys.monomers)
num_species(sys::FieldSystem) = length(sys.species)
volume(sys::FieldSystem) = sys.cell.volume

has_monomer(sys::FieldSystem, mid) = haskey(sys.monomers, mid)

function monomer_fractions!(sys::FieldSystem)
	bulk = sys.density_bulk
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
	add_monomer!(sys, mon)

Add a `Monomer` to the `FieldSystem`.
"""
function add_monomer!(sys::FieldSystem, mon::Monomer)
	if has_monomer(sys, mon)
		@warn "Replacing monomer with id = $(mon.id) in FieldSystem."
	end

	# Add to system, add new fields
	sys.monomers[mon.id] = mon
	sys.fields[mon.id] = zeros(Float64, sys.npw)
	sys.potentials[mon.id] = zeros(Float64, sys.npw)
	sys.residuals[mon.id] = zeros(Float64, sys.npw)
	sys.density[mon.id] = zeros(Float64, sys.npw)
	sys.density_bulk[mon.id] = 0.0

	return nothing
end

"""
	add_interaction!(sys, itx)

Add an `AbstractInteraction` to the `FieldSystem`.
"""
function add_interaction!(sys::FieldSystem, itx::AbstractInteraction)
	# Check to make sure has valid monomers
	for mid in itx.mids
        if !has_monomer(sys, mid)
            error("Interaction contains monomer with id = $(mid) not present in system.")
        end
    end

	setup!(itx, sys)
	push!(sys.interactions, itx)
	return nothing
end

"""
	add_species!(sys, species, phi/mu)

Add an `AbstractSpecies` to the `FieldSystem`.
If the system ensemble attribute is `Canonical`,
the third argument is interpreated as the bulk volume fraction of the species.
Conversely, this argument is the chemical potential of the species
for a `Grand` ensemble system.

Before relaxing the system configuration
the bulk volume fractions of all species in the system must sum to unity.
This constrains the possible values of `phi` and `mu`,
and is enforced in the program.
"""
function add_species!(sys::FieldSystem, species::AbstractSpecies, val::Real)
	# Check to make sure species has valid monomers
	for mid in species.mids
		if !has_monomer(sys, mid)
			error("Species contains monomer with id = $(mid) not present in system.")
		end
	end

	setup!(species, sys)
	push!(sys.species, species)
	# Add value based on ensemble of simulation
	if sys.ensemble == Canonical
		push!(sys.phi_species, val)
		push!(sys.mu_species, 0.0)
	elseif sys.ensemble == Grand
		push!(sys.phi_species, 0.0)
		push!(sys.mu_species, val)
	end

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
		error("System validation error: Fields must be initialized before SCFT simulation.")
	end

	# Check bulk volume fraction constraint
	if sys.ensemble == Canonical
		if !(sum(sys.phi_species) ≈ 1.0)
			error("System validation error: Species volume fractions do not sum to unity.")
		end
	end

	# Check bulk charge neutrality
	bulk_charge = 0.0
	monomer_fractions!(sys)
	for (mid, mon) in sys.monomers
		bulk_charge += sys.density_bulk[mid] * mon.charge / mon.size		
	end
	if !(bulk_charge ≈ 0.0)
		error("System validation error: Bulk system not charge neutral.")
	end

	# Check if any monomer types do not have associated species
	for mid in keys(sys.monomers)
		has_species = false
		for species in sys.species
			if has_monomer(species, mid)
				has_species = true
				break
			end
		end
		if !has_species
			error("System validation error: Monomer with id = $(mid) not found in any species.")
		end
	end

	return nothing
end

#==============================================================================#
# Methods for system density, potentials, and compressibility
#==============================================================================#

"""
	update_mu_phi!(sys)

Update species volume fractions and chemical potentials.
"""
function update_mu_phi!(sys::FieldSystem)
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
	update_density!(sys)

Compute density fields for all monomer types from all species in the system.
"""
function update_density!(sys::FieldSystem)
	# Zero density fields
	for rho in values(sys.density); rho .= 0.0; end

	# Loop over all species, update density on block/type level
	for species in sys.species
		update_density!(species)
	end

	# Update all bulk fractions/chemical potentials
	update_mu_phi!(sys)

	# Pool the updated density by type at the system level
	for (isp, species) in enumerate(sys.species)
		phi = sys.phi_species[isp]
		for mid in species.mids
			@. sys.density[mid] += phi * species.density[mid]
		end
	end

	# Sum all current density components for each grid point
	sys.density_sum .= 0.0
	for rho in values(sys.density)
		@simd for i in eachindex(rho)
			@inbounds sys.density_sum[i] += rho[i]
		end
	end

	return nothing
end

"""
	update_potentials!(sys)

Compute the differentials of the excess Hamiltonian w.r.t. density fields,
 ``\\frac{\\delta H_{exc}}{\\delta \\phi_\\alpha}``.
"""
function update_potentials!(sys::FieldSystem)
	# Zero current potentials
	for pot in values(sys.potentials)
		pot .= 0.0
	end

	# Add excess interactions
	for (mid, pot) in sys.potentials
		for itx in sys.interactions
			add_potential!(itx, mid, pot)
		end
	end

	return nothing
end 

"""
	update_eta!(sys)

Compute the Lagrange field ``\\eta`` based on the target scheme for enforcing incompressibility.
Assumes the system density fields and potentials are updated prior to calling.
"""
function update_eta!(sys::FieldSystem)
	comp = sys.comp

	# Compressible system
	if comp.inv_zeta > 0.0
		@. sys.eta = (1.0 / comp.inv_zeta) * (sys.density_sum - 1.0)
	end

	# Incompressible system
	if comp.incompressible
		sys.eta .= 0.0
		for (mid, omega) in sys.fields
			@. sys.eta += omega - sys.potentials[mid]
		end
		sys.eta ./= num_monomers(sys)
	end

	return nothing
end

"""
	update_stress!(sys)

Compute the free energy variation of the system for each cell parameter.
"""
function update_stress!(sys::FieldSystem)
	# Reset to zero
	sys.stress .= 0.0

	# Add stress from each contributing species
	for (isp, species) in enumerate(sys.species)
		phi = sys.phi_species[isp]
		sys.stress .-= phi * scf_stress(species)
	end

	return nothing
end

"""
	set_omega_uniform!(sys)

Sets the value of all chemical potential fields equal
to the summation of interaction potentials, 
corresponding to the mean-field limit with vanishing constraint fields.
"""
function set_omega_uniform!(sys::FieldSystem)
	update_potentials!(sys)
	for (mid, omega) in sys.fields
		@. omega = sys.potentials[mid]
		if sys.ensemble == Canonical; 
			omega .-= mean(omega); 
		end
	end
	return nothing	
end

"""
    make_residuals!(sys)

Compute the gradient of the system Hamiltonian w.r.t. density fields,
equal to the residual values for the system fields at the mean-field condition.
"""
function make_residuals!(sys::FieldSystem)
	update_density!(sys)
	update_potentials!(sys)
	update_eta!(sys)

    # Make residuals from fields and calculated potentials/eta field
	for (mid, res) in sys.residuals
		@. res = sys.potentials[mid] + sys.eta - sys.fields[mid]
		if sys.comp.incompressible
			@. res += sys.density_sum - 1.0
		end
        if sys.ensemble == Canonical; res .-= mean(res); end
	end
    
	return nothing
end

#==============================================================================#

"""
	free_energy(sys)

Return the intensive Helmholtz free energy for the interacting, inhomogeneous system.
"""
function free_energy(sys::FieldSystem)
	# Add translational entropy terms
	ftrans = 0.0
	for (isp, species) in enumerate(sys.species)
		phi = sys.phi_species[isp]
		mu = sys.mu_species[isp]
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
	ffield /= num_grid(sys)

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

Return the intensive Helmholtz free energy for a system at the specified bulk volume fractions.
"""
function free_energy_bulk(sys::FieldSystem)
	fhelm = 0.0

	# Update bulk fractions
	monomer_fractions!(sys)

	# Add translational entropy terms
	for (isp, species) in enumerate(sys.species)
		phi = sys.phi_species[isp]
		if phi > 1e-8
			fhelm += phi * (log(phi) - 1.0) / species.Nref
		end
	end

	# Add the excess interaction energy
	for itx in sys.interactions
		fhelm += energy_bulk(itx)
	end

	return fhelm
end