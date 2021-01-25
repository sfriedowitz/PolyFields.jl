"""
	EulerUpdater(; method =  :PECE, lam = 0.05)

A field updater implementing predictor-corrector Euler time stepping with static step size `lam`.
Available options for `method` include:

* PEC - predict-evaluate-correct scheme, with 1 density evaluation per step
* PECE - predict-evaluate-correct-evaluate scheme, with 2 density evaluations per step
"""
struct EulerUpdater <: AbstractFieldUpdater
	method :: Symbol
	lam    :: Float64
	temp   :: Dict{Int,FieldGrid{Float64}} # Used as pre-allocated storage for half-steps

	function EulerUpdater(; method::Symbol = :PECE, lam::Real = 0.05)
		if !(method in (:PEC, :PECE))
			error("Invalid `method` parameter for EulerUpdater. Options include: (:PEC, :PECE)")
		end
		return new(method, lam, Dict())
	end
end

#==============================================================================#
# Methods
#==============================================================================#

Base.show(io::IO, updater::EulerUpdater) = @printf(io, "EulerUpdater(method = %s, lam = %.3g)", updater.method, updater.lam)

function setup!(updater::EulerUpdater, sys::FieldSystem)
	copydict!(updater.temp, sys.fields)
	return nothing
end

#==============================================================================#

function update_state!(sys::FieldSystem, updater::EulerUpdater)
    # f(w_t): dH/dphi = dH_exc/dphi + eta - w
	# Time-stepping scheme:
	#	Predictor: wp_t+1 = w_t + dt * f(w_t)
	#	Corrector: w_t+1  = w_t + 0.5 * dt * (f(w_t) + f(wp_t+1))

	# 1) Predict
	for (mid, omega) in sys.fields
		@. omega += updater.lam * sys.residuals[mid]
		@. updater.temp[mid] = omega + 0.5 * updater.lam * sys.residuals[mid]
	end

	# 2 Evaluate
	residuals!(sys)

	# 3) Correct
	for (mid, omega_temp) in updater.temp
		@. omega_temp += 0.5 * updater.lam * sys.residuals[mid]
	end
	copydict!(sys.fields, updater.temp)

	# 4) Evaluate (PECE mode)
	if updater.method == :PECE
		residuals!(sys)
	end

	return nothing
end