"""
	set_interaction!(itx, args...)

Set interaction parameters for the whole interaction or specific monomers.
"""
set_interaction!(its::AbstractInteraction, args...) = nothing

"""
    energy(itx)

Return the energy of the `AbstractInteraction` due to the current system configuration.
"""
energy(itx::AbstractInteraction) = nothing

"""
    bulk_energy(itx)

Return the energy of the `AbstractInteraction` due to the homogeneous bulk composition.
"""
bulk_energy(itx::AbstractInteraction) = nothing

"""
    potential!(itx, alpha, potential)

Add the potential contribution for monomer ID `alpha` to the grid.
"""
potential!(itx::AbstractInteraction, alpha::Integer, pot::FieldGrid) = nothing