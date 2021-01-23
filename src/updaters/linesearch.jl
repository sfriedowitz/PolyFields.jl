#==============================================================================#
# Wrapped line-search routine
#==============================================================================#

"""
	perform_line_search(ls, a0, sys, steps, temp)

Perform a line-search provided by `ls` to determine the optimal step size
to reduce the system free energy in the direction of `steps`.
"""
function perform_line_search(ls, sys::FieldSystem, a0::Real, steps::Dict{TK,TV}, temp::Dict{TK,TV}) where {TK<:Integer,TV<:FieldGrid}
	# Univariate functions
	phi(a) = phi_ls(a, sys, steps, temp)
	dphi(a) = dphi_ls(a, sys, steps, temp)
	phi_dphi(a) = phi_dphi_ls(a, sys, steps, temp)

	fx, d0 = phi_dphi(a0)
	a, fx = ls(phi, dphi, phi_dphi, a0, fx, d0)

	return a
end

#==============================================================================#
# Single variable line search methods
#==============================================================================#

# IMPORTANT NOTES:
# 	The argument `steps` must not be a reference to system field arrays (i.e. residuals)
#	This is because we want to leave the system state unmodified for each trial evaluaiton of `a`
#	We thus allow the system density and resiudals to live in a modified state between evals at the new step point
#	This means that `steps` should not be a reference to `sys.residuals`, as it would be modified in the linesearch
	
function phi_ls(a::Real, sys::FieldSystem, steps::Dict{TK,TV}, temp::Dict{TK,TV}) where {TK<:Integer,TV<:FieldGrid}
	# Evaluate the objective for a given step
	# Steps:

	# 1) Copy system fields to storage
	copy_dict!(sys.fields, temp)

	# 2) Increment system fields by step
	for (mid, omega) in sys.fields
		@. omega += a * steps[mid]
	end

	# 3) Update system density with those fields
	make_residuals!(sys)

	# 3) Compute objective with new density profile
	#fx = free_energy(sys)
	fx = field_error(sys)

	# 4) Return system fields to original field state
	copy_dict!(temp, sys.fields)

	return fx
end

function dphi_ls(a::Real, sys::FieldSystem, steps::Dict{TK,TV}, temp::Dict{TK,TV}) where {TK<:Integer,TV<:FieldGrid}
	# Evaluate the gradient at the new location
	# Return the dot product of the gradient and the step
	# Steps:

	# 1) Copy system fields to storage
	copy_dict!(sys.fields, temp)

	# 2) Increment system fields by step
	for (mid, omega) in sys.fields
		@. omega += a * steps[mid]
	end

	# 2) Update system density with those fields
	make_residuals!(sys)

	# 3) Compute directional derivative dot(gradients, steps)
	d = dot_dicts(sys.residuals, steps)

	# 4) Return system fields to original field state
	copy_dict!(temp, sys.fields)

	return d
end

function phi_dphi_ls(a::Real, sys::FieldSystem, steps::Dict{TK,TV}, temp::Dict{TK,TV}) where {TK<:Integer,TV<:FieldGrid}
	# Return both objective and gradient dot product at new step
	# Steps:

	# 1) Copy fields to storage
	copy_dict!(sys.fields, temp)

	# 2) Increment system fields by step
	for (mid, omega) in sys.fields
		@. omega += a * steps[mid]
	end

	# 2) Update system density with those fields
	make_residuals!(sys)

	# 3) Compute objective at new location
	#fx = free_energy(sys)
	fx = field_error(sys)

	# 4) Compute directional derivative dot(residuals, steps)
	d = dot_dicts(sys.residuals, steps)

	# 5) Return system to original field state
	copy_dict!(temp, sys.fields)

	return (fx, d)
end