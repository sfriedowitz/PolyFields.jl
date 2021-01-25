"""
    mutable struct DomainUpdater <: AbstractFieldUpdater

    DomainUpdater(; nper = 5, nskip = 100, lam = 0.1, pmax = 5.0, rlxn = 0.9)

"""
mutable struct DomainUpdater <: AbstractFieldUpdater
    nper      :: Int
    nskip     :: Int
    lam       :: Float64
    pmax      :: Float64
    rlxn      :: Float64

    oldparams :: Vector{Float64}
    step      :: Vector{Float64} # Newton step
    stress    :: Vector{Float64} # First derivative of free energy
    dstress   :: Matrix{Float64} # Jacobian matrix for cell parameters
    system    :: Option{FieldSystem}

    function DomainUpdater(; nper::Integer = 10, nskip::Integer = 200,
        lam::Real = 0.1, pmax::Real = 1.0, rlxn::Real = 0.5)
        return new(nper, nskip, lam, pmax, rlxn, [], [], [], zeros(0,0), nothing)
    end
end

#==============================================================================#

Base.show(io::IO, domain::DomainUpdater) = @printf(io, "DomainUpdater(nper = %d, nskip = %d)", domain.nper, domain.nskip)

stress_error(domain::DomainUpdater) = norm(domain.stress)

function setup!(domain::DomainUpdater, sys::FieldSystem)
    domain.system = sys

    npr = nparams(sys.cell)
    domain.stress = scfstress(sys)
    domain.oldparams = zeros(npr)
    domain.step = zeros(npr)
    domain.dstress = zeros(npr, npr)

    return nothing
end

function scfstress!(domain::DomainUpdater)
    @assert !isnothing(domain.system)
    domain.stress .= scfstress(domain.system)
    return nothing
end

function step!(domain::DomainUpdater; nsteps::Integer = 1, tol::Real = 1e-5)
    @assert !isnothing(domain.system)
    sys = domain.system
    cell = sys.cell

    # Calculate current stress
    scfstress!(domain)
    err = stress_error(domain)

    converged = false
    iter = 0

    while iter < nsteps && !converged
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

        # Update to new stress
        update!(cell)
        density!(sys)
        residuals!(sys)
        scfstress!(domain)

        err = stress_error(domain)
        iter += 1
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