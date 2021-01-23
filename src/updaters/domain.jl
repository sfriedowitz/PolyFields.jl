"""
    DomainUpdater(; lambda = 1.0, skip = 50)

"""
mutable struct DomainUpdater
    tol         :: Float64
    skip        :: Int

    old_params  :: Vector{Float64}
    step        :: Vector{Float64} # Newton step
    stress      :: Vector{Float64} # First derivative of free energy
    dstress     :: Matrix{Float64} # Jacobian matrix for cell parameters
    system      :: Option{FieldSystem}

    function DomainUpdater(; tol::Real = 1e-4, skip::Integer = 1)
        return new(tol, skip, [], [], [], zeros(0,0), nothing)
    end
end

#==============================================================================#
# Constructors
#==============================================================================#

function DomainUpdater(sys::FieldSystem)
    nparams = num_cell_params(sys.cell)
    old_params = zeros(nparams)

    # Initial updater fields
    stress = scf_stress(sys)
    step = zeros(nparams)
    dstress = zeros(nparams, nparams)

    return DomainUpdater(sys, old_params, step, stress, dstress)
end

#==============================================================================#
# Methods
#==============================================================================#

Base.show(io::IO, domain::DomainUpdater) = @printf(io, "DomainUpdater(tol = %.3e, skip = %d)", domain.tol, domain.skip)

function setup!(domain::DomainUpdater, sys::FieldSystem)
    domain.system = sys
    nparams = num_cell_params(sys.cell)

    domain.stress = scf_stress(sys)
    domain.old_params = zeros(nparams)
    domain.step = zeros(nparams)
    domain.dstress = zeros(nparams, nparams)

    return nothing
end

stress_error(domain::DomainUpdater) = norm(domain.stress) 

function step_ready(domain::DomainUpdater, step::Integer, err::Real)
    return step % domain.skip == 0 && err < domain.tol
end

function relax_cell!(domain::DomainUpdater)
    sys = domain.system
    cell = sys.cell

    # Calculate stress, Jacobian, and step direction
    domain.stress .= scf_stress(sys)
    iter = 0

    while iter < 1
        stress_response!(domain)
        if cell.dim == 1
            domain.step .= -1.0 * domain.stress[1] / domain.dstress[1]
        else
            mul!(domain.step, inv(domain.dstress), domain.stress)
            domain.step .*= -1.0
        end

        # Rescale step-size if needed
        step_max = 0.5 * max(norm(cell.params), num_cell_params(cell))
        step_size = norm(domain.step)
        if step_size > step_max
            domain.step .*= step_max / step_size
        end

        # Step cell parameters, update cell and stress
        cell.params .+= domain.step
        update_cell!(cell)
        update_density!(sys)
        update_residuals!(sys)
        domain.stress .= scf_stress(sys)

        #println(domain.stress)

        iter += 1
    end

    return nothing
end

function stress_response!(domain::DomainUpdater)
    # Finite difference step size
    dx = 1e-8
    x = 1.0 / dx

    # Offload fields, save old parameters of cell
    sys = domain.system
    cell = sys.cell
    domain.old_params .= cell.params

    # For each cell param, calculate response to small perturbation
    for k = 1:num_cell_params(cell)
        cell.params .= domain.old_params
        cell.params[k] += dx
        update_cell!(cell)
        update_density!(sys)

        new_stress = scf_stress(sys)
        @. domain.dstress[:,k] = x * (new_stress - domain.stress)
    end

    # Return to original state
    cell.params .= domain.old_params
    update_cell!(cell)
    update_density!(sys)

    return nothing
end