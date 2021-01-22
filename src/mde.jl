"""
	struct PseudoSpectralSolver

An iteration scheme used for solving modified diffusion equations using pseudo-spectral Fourier methods.

Parameter `method` specifies whether to use a second-order (RK2) or fourth-order (RQM4) time-stepping scheme.

### Constructors
```julia
PseudoSpectralSolver(npw; method)
PseudoSpectralSolver(Nx, Ny, Nz; method)
```
"""
struct PseudoSpectralSolver
	npw    :: NTuple{3,Int}
	method :: String

	# Operator grids for solving MDE
	LW1    :: PWGrid{Float64}
	LW2    :: PWGrid{Float64}
	LD1    :: PWGrid{Float64}
	LD2    :: PWGrid{Float64}
end

#==============================================================================#
# Constructors
#==============================================================================#

function PseudoSpectralSolver(npw::NTuple{3,<:Integer}; method::AbstractString = "RK2")
	method = strip(uppercase(method))
	if method != "RK2" && method != "RQM4"
		throw(ArgumentError("Invalid update method `$(method)` for PseudoSpectralSolver."))
	end

	# Operator grids
	LW1 = zeros(Float64, npw)
	LW2 = zeros(Float64, npw)
	LD1 = zeros(Complex{Float64}, floor(Int, npw[1]/2 + 1), npw[2], npw[3])
	LD2 = zeros(Complex{Float64}, floor(Int, npw[1]/2 + 1), npw[2], npw[3])

	return PseudoSpectralSolver(npw, method, LW1, LW2, LD1, LD2)
end

PseudoSpectralSolver(Nx::Integer, Ny::Integer, Nz::Integer; kwargs...) = PseudoSpectralSolver((Nx, Ny, Nz); kwargs...)

#==============================================================================#
# Methods
#==============================================================================#

show(io::IO, solver::PseudoSpectralSolver) = @printf(io, "PseudoSpectralSolver(npw = %s)", solver.npw)

#==============================================================================#

"""
	update_operators!(solver, omega, ksq, size, b, ds)

Update differential operators for a given set of `omega` and `ksq` grids
with specified monomer size and chain length parameters.
"""
function update_operators!(solver::PseudoSpectralSolver, omega::PWGrid, ksq::PWGrid, vbar::Real, b::Real, ds::Real)
	@assert size(omega) == solver.npw
	@assert size(ksq) == size(solver.LD1)

	lw_coeff = ds / 2.0
	ld_coeff = ds * b^2 / 6.0

	@. solver.LW1 = exp(-lw_coeff * vbar * omega)
	@. solver.LW2 = exp(-lw_coeff * vbar * omega / 2.0)

	@. solver.LD1 = exp(-ld_coeff * ksq)
	@. solver.LD2 = exp(-ld_coeff * ksq / 2.0)

	return nothing
end

"""
	propagate!(solver, fft, qin, qout)

Take one step forward along the chain propagator from input `qin` to output `qout`.
Utilizes pre-allocated fast-Fourier transforms of the correct dimension from the provided `FFTBuddy`.
"""
function propagate!(solver::PseudoSpectralSolver, fft::FFTBuddy, qin::PWGrid, qout::PWGrid)
	@assert solver.npw == fft.npw

	# Offload params
	FT, iFT = fft.FT, fft.iFT
	LW1, LW2, LD1, LD2 = solver.LW1, solver.LW2, solver.LD1, solver.LD2

	# Pre-allocated temp grids from FFT
	q1, q2, qr = rgrid(fft, 1), rgrid(fft, 2), rgrid(fft, 3)
	qk = kgrid(fft, 1)

	if solver.method == "RK2"
		@. qr = LW1 * qin
		mul!(qk, FT, qr)
		@. qk *= LD1
		mul!(qr, iFT, qk)
		@. qout = LW1 * qr

	elseif solver.method == "RQM4"
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