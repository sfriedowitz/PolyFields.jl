"""
	set_interaction(itx, args...)

Set coefficients for a specific interaction between chemical species.
"""
set_interaction(its::AbstractInteraction, args...) = nothing

"""
    energy(itx)

Return the energy contribution of the `AbstractInteraction` for a given set of system fields.
"""
energy(itx::AbstractInteraction) = nothing

"""
    bulk_energy(itx)

Return the energy contribution of the interaction corresponding 
to the homogeneous bulk composition.
"""
bulk_energy(itx::AbstractInteraction) = nothing

"""
    potential!(itx, alpha, grid)

Add the potential contribution for monomer `alpha` to the provided field grid.
"""
potential!(itx::AbstractInteraction, alpha::Integer, grid::FieldGrid) = nothing