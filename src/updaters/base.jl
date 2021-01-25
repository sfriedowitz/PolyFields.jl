"""
    step!(updater)

Perform a single step to update the system fields 
using the methods defined by a provided `AbstractFieldUpdater`.
"""
step!(updater::AbstractFieldUpdater) = nothing

setup!(species::AbstractFieldUpdater, sys::FieldSystem) = nothing