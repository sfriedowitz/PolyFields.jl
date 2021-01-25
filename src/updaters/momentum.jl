"""
	mutable struct MomentumUpdater <: AbstractFieldUpdater

	MomentumUpdater(; method = :MOM, lam = 0.01, gamma = 0.9, beta1 = 0.9, beta2 = 0.999)

A field updater implementing various momentum-based iteration techniques.
These techniques track past gradient and step directions
to obtain a better current step direction.

Options for method include:
* MOM - Standard momentum update
* RMS - RMSprop update
* ADAM - Combination of ADAgrad and RMSprop, used commonly in NN training
* NGS - Nesterov accelerated gradient
"""
mutable struct MomentumUpdater <: AbstractFieldUpdater
	method :: Symbol
	lam    :: Float64
	gamma  :: Float64
	beta1  :: Float64
	beta2  :: Float64

	step   :: Dict{Int,FieldGrid{Float64}} # Final step for each method
	mt     :: Dict{Int,FieldGrid{Float64}} # First moment of gradient moving average
	vt     :: Dict{Int,FieldGrid{Float64}} # Second moment of gradient moving average
	temp   :: Dict{Int,FieldGrid{Float64}}

	system :: Option{FieldSystem}

	function MomentumUpdater(; method::Symbol = :MOM, lam::Real = 0.01, 
		gamma::Real = 0.9, beta1::Real = 0.9, beta2::Real = 0.999)
		@assert 0.0 < gamma < 1.0
		@assert 0.0 < beta1 < 1.0
		@assert 0.0 < beta2 < 1.0

		if !(method in (:MOM, :RMS, :ADAM, :NAG))
			error("Invalid `method` parameter for MomentumUpdater. Options include: (:MOM, :RMS, :ADAM, :NAG)")
		end
		return new(method, lam, gamma, beta1, beta2, Dict(), Dict(), Dict(), Dict(), nothing)
	end
end

#==============================================================================#

Base.show(io::IO, updater::MomentumUpdater) = @printf(io, "MomentumUpdater(method = %s, lam = %.3g)", updater.method, updater.lam)

function setup!(updater::MomentumUpdater, sys::FieldSystem)
	updater.system = sys

	# Copy fields to step and temp
	copydict!(updater.step, sys.residuals)
	copydict!(updater.temp, sys.fields)

	# Initialize first and second moments of gradient
	for mid in keys(sys.fields)
		updater.mt[mid] = zeros(Float64, sys.dims)
		updater.vt[mid] = zeros(Float64, sys.dims)
	end

	return nothing
end

function step!(updater::MomentumUpdater)
	@assert !isnothing(updater.system)
	sys = updater.system

	# Update steps
	if updater.method == :MOM
		_step_mom!(updater)
	elseif updater.method == :RMS
		_step_rms!(updater)
	elseif updater.method == :ADAM
		_step_adam!(updater)
	elseif updater.method == :NAG
		_step_nag!(updater)
	end

	# Take step
	for (mid, omega) in sys.fields
		@. omega += updater.lam * updater.step[mid]
	end

	return nothing
end

function _step_mom!(updater::MomentumUpdater)
	sys = updater.system
	residuals!(sys)

	# Update the first moment history
	for (mid, m) in updater.mt
		@. m = updater.beta1 * m + (1 - updater.beta1) * sys.residuals[mid]
	end

	# Copy updated momentum to step
	copydict!(updater.step, updater.mt)

	return nothing
end

function _step_rms!(updater::MomentumUpdater)
	sys = updater.system
	residuals!(sys)

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

function _step_adam!(updater::MomentumUpdater)
	sys = updater.system
	residuals!(sys)

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

function _step_nag!(updater::MomentumUpdater)
	sys = updater.system

	# Store current fields in temp
	copydict!(updater.temp, sys.fields)

	# Take a step with previous momentum term
	for (mid, omega) in sys.fields
		@. omega += updater.lam * updater.mt[mid]
	end

	# Compute gradient at update location
	residuals!(sys)

	# Update momentum w/ new residuals
	for (mid, m) in updater.mt
		@. m = updater.gamma * m + (1 - updater.gamma) * sys.residuals[mid]
	end

	# Copy updated momentum to step
	copydict!(updater.step, updater.mt)

	# Revert fields back to their previous values
	copydict!(sys.fields, updater.temp)

	return nothing
end