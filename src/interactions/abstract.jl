#==============================================================================#
# Methods for abstract interactions
#==============================================================================#

"""
    set_interaction!(itx, args...)

Set interaction parameters specific to a given `AbstractInteraction`.
"""
set_interaction!(itx::AbstractInteraction, args...) = nothing

"""
    energy(itx)

Return the energy contribution of the `AbstractInteraction` for a given set of system fields.
"""
energy(itx::AbstractInteraction) = nothing

"""
    energy_bulk(itx)

Return the energy contribution of the `AbstractInteraction` corresponding 
to the homogeneous bulk composition.
"""
energy_bulk(itx::AbstractInteraction) = nothing

"""
    add_potential!(itx, alpha, pot)

Add the potential contribution for species `alpha` to the array grid `pot`.
"""
add_potential!(itx::AbstractInteraction, alpha::Integer, pot::NPWGrid) = nothing