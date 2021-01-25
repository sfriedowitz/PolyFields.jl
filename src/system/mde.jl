"""
	struct PseudoSpectralSolver <: AbstractPropagator

An iteration scheme used for solving modified diffusion equations using pseudo-spectral Fourier methods.
The `method` parameter specifies whether to use a second-order (:RK2) or fourth-order (:RQM4) time-stepping scheme.
"""
struct PseudoSpectralSolver<: AbstractPropagator
	dims   :: NTuple{3,Int}
	method :: Symbol

	# Operator grids for solving MDE
	LW1    :: FieldGrid{Float64}
	LW2    :: FieldGrid{Float64}
	LD1    :: FieldGrid{Float64}
	LD2    :: FieldGrid{Float64}
end

function PseudoSpectralSolver(dims::NTuple{3,<:Integer}; method::Symbol = :RK2)
	if method != :RK2 && method != :RQM4
		throw(ArgumentError("Invalid propagator method `$(method)` for PseudoSpectralSolver."))
	end

	# Operator grids
	LW1 = zeros(Float64, dims)
	LW2 = zeros(Float64, dims)
	LD1 = zeros(Float64, ksize(dims))
	LD2 = zeros(Float64, ksize(dims))

	return PseudoSpectralSolver(dims, method, LW1, LW2, LD1, LD2)
end

PseudoSpectralSolver(Nx::Integer, Ny::Integer, Nz::Integer; kwargs...) = PseudoSpectralSolver((Nx, Ny, Nz); kwargs...)

#==============================================================================#

Base.show(io::IO, solver::PseudoSpectralSolver) = @printf(io, "PseudoSpectralSolver(dims = %s)", solver.dims)

"""
	update!(solver, omega, ksq, vol, b, ds)

Update differential operators for a given set of `omega` and `ksq` grids.
"""
function update!(solver::PseudoSpectralSolver, omega::FieldGrid, ksq::FieldGrid, vol::Real, b::Real, ds::Real)
	@assert size(omega) == solver.dims
	@assert size(ksq) == size(solver.LD1)

	lw_coeff = ds / 2.0
	ld_coeff = ds * b^2 / 6.0

	@. solver.LW1 = exp(-lw_coeff * vol * omega)
	@. solver.LW2 = exp(-lw_coeff * vol * omega / 2.0)

	@. solver.LD1 = exp(-ld_coeff * ksq)
	@. solver.LD2 = exp(-ld_coeff * ksq / 2.0)

	return nothing
end

"""
	propagate!(solver, fft, qin, qout)

Take one step forward along the chain propagator from input `qin` to output `qout`.
Utilizes pre-allocated fast-Fourier transforms of the correct dimension from the provided `FFTHolder`.
"""
function propagate!(solver::PseudoSpectralSolver, plan::FFTHolder, qin::FieldGrid, qout::FieldGrid)
	@assert solver.dims == plan.dims

	# Offload params
	FT, iFT = plan.FT, plan.iFT
	LW1, LW2, LD1, LD2 = solver.LW1, solver.LW2, solver.LD1, solver.LD2

	# Pre-allocated temp grids from FFT
	q1, q2, qr = rgrid(plan, 1), rgrid(plan, 2), rgrid(plan, 3)
	qk = kgrid(plan, 1)

	if solver.method == :RK2
		@. qr = LW1 * qin
		mul!(qk, FT, qr)
		@. qk *= LD1
		mul!(qr, iFT, qk)
		@. qout = LW1 * qr

	elseif solver.method == :RQM4
		# Full step
		@. qr = LW1 * qin
		mul!(qk, FT, qr)
		@. qk *= LD1
		mul!(qr, iFT, qk)
		@. q1 = LW1 * qr
		
		# Two half-steps
		@. qr = LW2 * qin
		mul!(qk, FT, qr)
		@. qk *= LD2
		mul!(qr, iFT, qk)
		@. qr *= LW1
		mul!(qk, FT, qr)
		@. qk *= LD2
		mul!(qr, iFT, qk)
		@. q2 = LW2 * qr

		# Richardson extrapolation
		@. qout = (4.0*q2 - q1)/3.0
	end

	return nothing
end