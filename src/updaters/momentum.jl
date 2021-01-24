"""
	MomentumUpdater(; method = "MOM", lam = 0.05, gamma = 0.9, beta1 = 0.9, beta2 = 0.999)


"""
mutable struct MomentumUpdater <: AbstractFieldUpdater
	method :: String
	lam    :: Float64
	gamma  :: Float64
	beta1  :: Float64
	beta2  :: Float64

	step   :: Dict{Int,FieldGrid{Float64}} # Final step for each method
	mt     :: Dict{Int,FieldGrid{Float64}} # First moment of gradient moving average
	vt     :: Dict{Int,FieldGrid{Float64}} # Second moment of gradient moving average
	temp   :: Dict{Int,FieldGrid{Float64}}

	function MomentumUpdater(;  
		method::AbstractString = "MOM", 
		lam::Real = 0.05, 
		gamma::Real = 0.9, 
		beta1::Real = 0.9, 
		beta2::Real = 0.999
		)
		@assert 0.0 < gamma < 1.0
		@assert 0.0 < beta1 < 1.0
		@assert 0.0 < beta2 < 1.0

		method = strip(uppercase(method))
		if !(method in ("MOM", "RMS", "ADAM", "NAG"))
			error("Invalid `method` parameter for MomentumUpdater. Options include: (MOM, RMS, ADAM, NAG)")
		end
		return new(method, lam, gamma, beta1, beta2, Dict(), Dict(), Dict(), Dict())
	end
end

#==============================================================================#
# Methods
#==============================================================================#

Base.show(io::IO, updater::MomentumUpdater) = @printf(io, "MomentumUpdater(method = %s, lam = %.3g)", updater.method, updater.lam)

function setup!(updater::MomentumUpdater, sys::FieldSystem)
	# Copy fields to step and temp
	copy_dict!(sys.residuals, updater.step)
	copy_dict!(sys.fields, updater.temp)

	# Initialize first moments of gradient
	for (mid, _) in sys.fields
		updater.mt[mid] = zeros(sys.dims)
	end

	# Initialize second moments of gradient
	for (mid, _) in sys.fields
		updater.vt[mid] = zeros(sys.dims)
	end

	return nothing
end

#==============================================================================#

function update_state!(sys::FieldSystem, updater::MomentumUpdater)
	# Update steps
	if updater.method == "MOM"
		_update_step_mom!(updater, sys)
	elseif updater.method == "RMS"
		_update_step_rms!(updater, sys)
	elseif updater.method == "ADAM"
		_update_step_adam!(updater, sys)
	elseif updater.method == "NAG"
		_update_step_nag!(updater, sys)
	end

	# Take step
	for (mid, omega) in sys.fields
		@. omega += updater.lam * updater.step[mid]
	end

	return nothing
end

function _update_step_mom!(updater::MomentumUpdater, sys::FieldSystem)
	# Compute new residuals
	make_residuals!(sys)

	# Update the first moment history
	for (mid, m) in updater.mt
		@. m = updater.beta1 * m + (1 - updater.beta1) * sys.residuals[mid]
	end

	# Copy updated momentum to step
	copy_dict!(updater.mt, updater.step)

	return nothing
end

function _update_step_rms!(updater::MomentumUpdater, sys::FieldSystem)
	# Compute new residuals
	make_residuals!(sys)

	# Update the second moment history
	for (mid, v) in updater.vt
		@. v = updater.beta2 * v + (1 - updater.beta2) * sys.residuals[mid]^2
	end

	# Update the step given new residuals and momentum
	for (mid, s) in updater.step
		@. s = sys.residuals[mid] / (sqrt(updater.vt[mid]) + 1e-8)
	end

	return nothing
end

function _update_step_adam!(updater::MomentumUpdater, sys::FieldSystem)
	# Compute new residuals
	make_residuals!(sys)

	# Update the first moment history
	for (mid, m) in updater.mt
		@. m = updater.beta1 * m + (1 - updater.beta1) * sys.residuals[mid]
	end

	# Update the second moment history
	for (mid, v) in updater.vt
		@. v = updater.beta2 * v + (1 - updater.beta2) * sys.residuals[mid]^2
	end

	# Compute the new step
	for (mid, s) in updater.step
		 @. s = updater.mt[mid] / (sqrt(updater.vt[mid]) + 1e-8)
	end

	return nothing
end

function _update_step_nag!(updater::MomentumUpdater, sys::FieldSystem)
	# Store current fields in temp
	copy_dict!(sys.fields, updater.temp)

	# Take a step with previous momentum term
	for (mid, omega) in sys.fields
		@. omega += updater.lam * updater.mt[mid]
	end

	# Compute gradient at update location
	make_residuals!(sys)

	# Update momentum w/ new residuals
	for (mid, m) in updater.mt
		@. m = updater.gamma * m + (1 - updater.gamma) * sys.residuals[mid]
	end

	# Copy updated momentum to step
	copy_dict!(updater.mt, updater.step)

	# Revert fields back to their previous values
	copy_dict!(updater.temp, sys.fields)

	return nothing
end