"""
    Compressibility(; inv_zeta = 0.0)

Construct an incompressibility constraint for the system. 
Parameter `inv_modulus` controls how incompressibility is treated.

* ``\\zeta^{-1} = 0.0``: System is treated as incompressible
* ``\\zeta^{-1} > 0.0``: System is treated as compressible with bulk modulus ``\\zeta``
* ``\\zeta^{-1} < 0.0``: No incompressibility is enforced
"""
struct Compressibility <: AbstractConstraint
    inv_zeta       :: Float64
    incompressible :: Bool
    field          :: FieldGrid{Float64}
    system         :: Option{FieldSystem}

    function Compressibility(; inv_modulus::Real = 0.0)
        incompressible = inv_modulus ≈ 0.0 ? true : false
        return new(inv_modulus, incompressible, [], nothing)
    end
end

"""
    Electroneutrality(; inv_modulus = 0.0)
"""
struct Electroneutrality <: AbstractConstraint
    inv_modulus    :: Float64
    electroneutral :: Bool
    field          :: FieldGrid{Float64}
    system         :: Option{FieldSystem}

    function Electroneutrality(; inv_modulus::Real = 0.0)
        electroneutral = inv_modulus ≈ 0.0 ? true : false
        return new(inv_modulus, electroneutral)
    end
end

#==============================================================================#

function setup!(c::AbstractConstraint, sys::FieldSystem)
    c.system = sys
    c.field = zeros(sys.npw)
    return nothing
end