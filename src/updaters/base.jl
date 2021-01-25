"""
    step!(updater)

Perform a single step to update the system fields 
using the methods defined by a provided `AbstractFieldUpdater`.
"""
step!(updater::AbstractFieldUpdater) = nothing

setup!(species::AbstractFieldUpdater, sys::FieldSystem) = nothing

"""
	updater_with_system(type, sys; kwargs...)

Return a field updater allocated for the provided system.
"""
function updater_with_system(TU::Type{<:AbstractFieldUpdater}, sys::FieldSystem; kwargs...)
	updater = TU(; kwargs...)
	setup!(updater, sys)
	return updater
end