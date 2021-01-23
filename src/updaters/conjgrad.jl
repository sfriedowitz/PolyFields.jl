"""
	ConjGradUpdater(; )
"""
mutable struct ConjGradUpdater <: AbstractFieldUpdater
	ls   :: Any
	nu   :: Float64
	lmin :: Float64
	lrat :: Float64

	lam  :: Float64
	d0   :: Float64

    pv   :: Dict{Int,FieldGrid{Float64}} # Current steps
	res0 :: Dict{Int,FieldGrid{Float64}} # Previous gradients
	temp :: Dict{Int,FieldGrid{Float64}} # Temp storage of fields during line-search

	function ConjGradUpdater(; ls, lam::Real = 0.1, nu::Real = 0.1, lmin::Real = 1e-3, lrat::Real = 10.0)
	 	return new(ls, nu, lmin, lrat, lam, 0.0, Dict(), Dict(), Dict())
	end
end

#==============================================================================#
# Methods
#==============================================================================#

Base.show(io::IO, updater::ConjGradUpdater) = @printf(io, "ConjGradUpdater(lam = %.3g)", updater.lam)

function setup!(updater::ConjGradUpdater, sys::FieldSystem)
	copy_dict!(sys.residuals, updater.pv)
	copy_dict!(sys.residuals, updater.res0)
	copy_dict!(sys.fields, updater.temp)

	# Calculate function value and directional derivative w/ initial fields
	updater.d0 = dot_dicts(updater.pv, sys.residuals)
	#updater.lam = max(0.01, 1.0 / (1 - updater.d0))

	return nothing
end

#==============================================================================#

function pr_update(curr::Dict{TK,TV}, prev::Dict{TK,TV}, epsilon = eps(Float64)) where {TK<:Integer,TV<:FieldGrid}
	num = dot_dicts(curr, curr) - dot_dicts(curr, prev)
	denom = dot_dicts(prev, prev) + epsilon
	return num / denom
end

function update_state!(sys::FieldSystem, updater::ConjGradUpdater)
	# Line search with current step directions
	lam_ls = updater.lam
	try
		lam_ls = perform_line_search(updater.ls, sys, 0.5, updater.pv, updater.temp)
		lam_ls = max(lam_ls, updater.lmin)
	catch err
		@warn "Error performing line-search: $(err)"
	end

	println(lam_ls)


	# Update fields with new step size
	for (mid, omega) in sys.fields
		@. omega += lam_ls * updater.pv[mid]
	end

	# Compute the new density and residuals
	make_residuals!(sys)

	# Determine new step direction: h_j = beta * h_j-1 - grad(F)
	# beta determined using Polak-Ribiere method
	beta = pr_update(sys.residuals, updater.res0)
	for (mid, step) in updater.pv
		@. step = beta*step + sys.residuals[mid]
	end

	# Test for orthogonality of the new and old search direction
	# Dot product of current (after step) and previous gradients
	ortho_val = abs(dot_dicts(sys.residuals, updater.res0) / dot_dicts(sys.residuals, sys.residuals))

	# Backup old values
	d_old = updater.d0
	copy_dict!(sys.residuals, updater.res0)
	updater.d0 = dot_dicts(updater.pv, sys.residuals)

	# # Check if CG reset to SD direction
	# reset_cg = false
	# if updater.d0 > 0.0 # Test for descent direction of line search
	# 	@warn "New CG direction is not a descent direction -- resetting to steepest descent."
	# 	reset_cg = true
	# elseif ortho_val > updater.nu # Test for sufficient orthogonality
	# 	@warn "New CG direction is not sufficiently orthogonal -- resetting to steepest descent."
	# 	reset_cg = true
	# elseif beta < 0.0 # Test for negative PR method
	# 	@warn "Poblak-Ribiere `beta` is negative -- resetting to steepest descent."
	# 	reset_cg = true
	# end

	# Reset search direction to gradient if necessary
	if beta < 0.0
		copy_dict!(sys.residuals, updater.pv)
		updater.d0 = dot_dicts(updater.pv, sys.residuals)
	end

	# # Determine new initial step size
	# slope_ratio = d_old / (updater.d0 - eps(Float64))
	# updater.lam = updater.lam * min(updater.lrat, slope_ratio)

	return nothing
end
