"""
    step!(updater)

Perform a single step to update the system fields 
using the methods defined by a provided `AbstractFieldUpdater`.
"""
step!(updater::AbstractFieldUpdater) = nothing

"""
    setup!(updater, sys)

Setup internal state of the field-updater for the provided system.
"""
setup!(species::AbstractFieldUpdater, sys::FieldSystem) = nothing