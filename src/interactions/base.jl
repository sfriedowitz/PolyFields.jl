"""
	set_interaction!(itx, args...)

Set interaction parameters for the whole interaction or specific monomers.
"""
set_interaction!(itx::AbstractInteraction, args...) = nothing
set_interaction!(itx::AbstractInteraction, mon::Monomer, args...) = set_interaction!(itx, mon.id, args...)
set_interaction!(itx::AbstractInteraction, mon1::Monomer, mon2::Monomer, args...) = set_interaction!(itx, mon1.id, mon2.id, args...)

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

setup!(itx::AbstractInteraction, sys::AbstractSystem)  = nothing
