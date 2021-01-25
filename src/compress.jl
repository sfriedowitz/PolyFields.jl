"""
    mutable struct Compressibility <: AbstractConstraint

    Compressibility(inv_zeta = 0.0)

A local incompressibility constraint on each grid point in a system.

Parameter `inv_zeta` controls how incompressibility is treated:
* ``ζ^{-1} = 0.0``: System is treated as incompressible
* ``ζ^{-1} > 0.0``: System is treated as compressible with bulk modulus ``ζ``
* ``ζ^{-1} < 0.0``: Incompressibility is not enforced at all
"""
mutable struct Compressibility <: AbstractConstraint
    inv_zeta :: Float64
    strict   :: Bool
    field    :: FieldGrid{Float64}
    rhosum   :: FieldGrid{Float64}
    system   :: Option{AbstractSystem}

    function Compressibility(inv_zeta::Real = 0.0)
        strict = inv_zeta ≈ 0.0 ? true : false
        return new(inv_zeta, strict, zeros(0,0,0), zeros(0,0,0), nothing)
    end
end

#==============================================================================#

Base.show(io::IO, c::Compressibility) = @printf(io, "Compressibility(inv_zeta = %.3f)", c.inv_zeta)

function setup!(c::Compressibility, sys::AbstractSystem)
    c.system = sys
    c.field = zeros(Float64, sys.dims)
    c.rhosum = zeros(Float64, sys.dims)
    return nothing
end

function update!(c::Compressibility)
    @assert !isnothing(c.system)
    sys = c.system

    # Update density sum
    c.rhosum .= 0.0
    for (mid, rho) in sys.density
        @. c.rhosum += rho
    end

    if c.strict
        # Incmpressible system
        @. c.field = c.rhosum - 1.0
        for (mid, omega) in sys.fields
            @. c.field += omega - sys.potentials[mid]
        end
        c.field ./= nmonomers(sys)
    elseif c.inv_zeta > 0.0
        # Compressible system
        @. c.field = (1.0/c.inv_zeta) * (c.rhosum - 1.0)
    end

    return nothing
end
