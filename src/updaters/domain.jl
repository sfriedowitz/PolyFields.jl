"""
    mutable struct DomainUpdater <: AbstractFieldUpdater

    DomainUpdater(; lambda = 1.0, skip = 50)

"""
mutable struct DomainUpdater <: AbstractFieldUpdater
    nsteps    :: Int
    nskip     :: Int
    verbose   :: Bool
    lam       :: Float64
    tol       :: Float64
    pmax      :: Float64
    rlxn      :: Float64

    oldparams :: Vector{Float64}
    step      :: Vector{Float64} # Newton step
    stress    :: Vector{Float64} # First derivative of free energy
    dstress   :: Matrix{Float64} # Jacobian matrix for cell parameters
    system    :: Option{FieldSystem}

    function DomainUpdater(; nsteps::Integer = 50, nskip::Integer = 1, verbose::Bool = true,
        lam::Real = 0.001, tol::Real = 1e-5, pmax::Real = 5.0, rlxn::Real = 0.9)
        return new(nsteps, nskip, verbose, lam, tol, pmax, rlxn, [], [], [], zeros(0,0), nothing)
    end
end

#==============================================================================#

Base.show(io::IO, domain::DomainUpdater) = @printf(io, "DomainUpdater(tol = %.1e, nskip = %d)", domain.tol, domain.nskip)

function setup!(domain::DomainUpdater, sys::FieldSystem)
    domain.system = sys

    npr = nparams(sys.cell)
    domain.stress = scfstress(sys)
    domain.oldparams = zeros(npr)
    domain.step = zeros(npr)
    domain.dstress = zeros(npr, npr)

    return nothing
end

function step!(domain::DomainUpdater)
    @assert !isnothing(domain.system)
    sys = domain.system
    cell = sys.cell

    # Calculate stress, Jacobian, and step direction
    domain.stress .= scfstress(sys)
    iter = 0
    err = norm(domain.stress)
    converged = err < domain.tol

    if converged
        return nothing
    elseif domain.verbose
        @printf("Relaxing cell for %d steps from initial stress norm: |σ| = %.3e.\n", iter, err)
    end

    while iter < domain.nsteps && !converged
        dstress!(domain)
        if cell.dim == 1
            domain.step .= -1.0 * domain.stress[1] / domain.dstress[1]
        else
            mul!(domain.step, inv(domain.dstress), domain.stress)
            domain.step .*= -1.0
        end

        # Rescale step-size if needed
        step_size = norm(domain.step)
        if step_size > domain.pmax
            domain.step .*= domain.pmax / step_size
        end

        # Step cell parameters, update cell and stress
        @. cell.params += domain.lam * domain.step
        @. cell.params = domain.rlxn*cell.params + (1 - domain.rlxn)*domain.oldparams

        update!(cell)
        density!(sys)
        residuals!(sys)

        domain.stress .= scfstress(sys)

        iter += 1
        err = norm(domain.stress)
        converged = err < domain.tol
    end

    if converged
        @printf("Converged on step %d with stress norm: |σ| = %.3e.\n", iter, err)
    else
        @printf("Finished relaxation of %d steps with stress norm: |σ| = %.3e.\n", iter, err)
    end

    return nothing
end

function dstress!(domain::DomainUpdater)
    @assert !isnothing(domain.system)
    sys = domain.system
    cell = sys.cell
    domain.oldparams .= cell.params

    # Finite difference step size
    dx = 1e-8
    x = 1.0 / dx

    # For each cell param, calculate response to small perturbation
    for k = 1:nparams(cell)
        cell.params .= domain.oldparams
        cell.params[k] += dx

        update!(cell)
        density!(sys)
        pstress = scfstress(sys)

        @. domain.dstress[:,k] = x * (pstress - domain.stress)
    end

    # Return to original state
    cell.params .= domain.oldparams
    density!(sys)

    return nothing
end