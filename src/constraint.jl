#==============================================================================#
# Abstract methods
#==============================================================================#

function setup!(c::AbstractConstraint, sys::FieldSystem)
    c.system = sys
    c.field = zeros(sys.npw)
    return nothing
end

#==============================================================================#
# Incompressibility
#==============================================================================#

"""
    Incompressibility(; inv_kappa = 0.0)

Construct an incompressibility constraint for the system. 
Parameter `inv_kappa` controls how incompressibility is treated.

* ``\\kappa^{-1} = 0.0``: System is treated as incompressible
* ``\\kappa^{-1} > 0.0``: System is treated as compressible with bulk modulus ``\\zeta``
* ``\\kappa^{-1} < 0.0``: No incompressibility is enforced
"""
struct Incompressibility <: AbstractConstraint
    inv_kappa      :: Float64
    incompressible :: Bool

    field          :: NPWGrid{Float64}
    system         :: Option{FieldSystem}

    function Compressibility(; inv_kappa::Real = 0.0)
        incompressible = inv_kappa ≈ 0.0 ? true : false
        return new(inv_kappa, incompressible, [], nothing)
    end
end

#==============================================================================#
# Electroneutrality
#==============================================================================#

"""
    Electroneutrality(; inv_kappa = 0.0)
"""
struct Electroneutrality <: AbstractConstraint
    inv_kappa      :: Float64
    electroneutral :: Bool

    field          :: NPWGrid{Float64}
    system         :: Option{FieldSystem}

    function Electroneutrality(; inv_kappa::Real = 0.0)
        electroneutral = kappa_inv ≈ 0.0 ? true : false
        return new(inv_kappa, electroneutral)
    end
end