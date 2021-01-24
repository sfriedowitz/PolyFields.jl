function setup!(cons::AbstractConstraint, sys::AbstractSystem)
    c.system = sys
    c.field = zeros(Float64, sys.dims)
    return nothing
end

#==============================================================================#


"""
    Compressibility(; inv_zeta = 0.0)

A local incompressibility constraint on each grid point in a system.

Parameter `inv_zeta` controls how incompressibility is treated:
* ``\\zeta^{-1} = 0.0``: System is treated as incompressible
* ``\\zeta^{-1} > 0.0``: System is treated as compressible with bulk modulus ``\\zeta``
* ``\\zeta^{-1} < 0.0``: No incompressibility is enforced
"""
mutable struct Compressibility <: AbstractConstraint
    inv_zeta :: Float64
    strict   :: Bool
    field    :: FieldGrid{Float64}

    function Compressibility(; inv_modulus::Real = 0.0)
        incompressible = inv_zeta ≈ 0.0 ? true : false
        return new(inv_zeta, incompressible, [])
    end
end

# """
#   update_eta!(sys)

# Compute the Lagrange field ``\\eta`` based on the target scheme for enforcing incompressibility.
# Assumes the system density fields and potentials are updated prior to calling.
# """
# function update_eta!(sys::FieldSystem)
#   comp = sys.comp

#   # Compressible system
#   if comp.inv_zeta > 0.0
#       @. sys.eta = (1.0 / comp.inv_zeta) * (sys.density_sum - 1.0)
#   end

#   # Incompressible system
#   if comp.incompressible
#       sys.eta .= 0.0
#       for (mid, omega) in sys.fields
#           @. sys.eta += omega - sys.potentials[mid]
#       end
#       sys.eta ./= num_monomers(sys)
#   end

#   return nothing
# end

#==============================================================================#

"""
    Electroneutrality(; inv_zeta = 0.0)
"""
mutable struct Electroneutrality <: AbstractConstraint
    inv_zeta :: Float64
    strict   :: Bool
    field    :: FieldGrid{Float64}

    function Electroneutrality(; inv_modulus::Real = 0.0)
        electroneutral = inv_zeta ≈ 0.0 ? true : false
        return new(inv_zeta, electroneutral)
    end
end