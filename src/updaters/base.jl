"""
    step!(sys, updater)

Perform a single step to update the system fields 
using the methods defined by a provided `AbstractFieldUpdater`.
"""
step!(updater::AbstractFieldUpdater) = nothing

"""
    field_error(sys)

Compute the error of the current residuals and the field configurations.
Follows the definition in the Matsen (2009) Anderson mixing papers:

``
\\epsilon = [\\sum_{\\alpha,i} r_{\\alpha,i}^2 / \\sum_{\\alpha,i} \\omega_{\\alpha,i}^2]^{1/2}
``
"""
function field_error(sys::FieldSystem)
    # Residual weights
    rsum = 0.0
    for (mid, res) in sys.residuals
        @simd for i in eachindex(res)
            @inbounds rsum += res[i]^2
        end
    end
    # Field weights
    wsum = eps(Float64)
    for (mid, omega) in sys.fields
        @simd for i in eachindex(omega)
            @inbounds wsum += omega[i]^2
        end
    end
    return sqrt(rsum) / sqrt(wsum)
end

"""
    stress_error(sys)

Compute the error related to stress in the system as the absolute value of the maximum stress component.
"""
function stress_error(sys::FieldSystem)
	max_val = 0.0
	for s in sys.stress
		if abs(s) > max_val
			max_val = abs(s)
		end
	end
	return max_val
end