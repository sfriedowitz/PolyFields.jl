"""
	struct FFTBuddy

A helper struct that stores pre-allocated utilities
for performing fast-Fourier transforms for a given grid-size.
Allocates `ngrids` temporary complex `NPWGrid`s for use as
temporary storage for in-place FFT operations.

### Constructors
```julia
FFTBuddy(npw; ngrids = 0, nthreads = 1)
```
"""
struct FFTBuddy
    npw    :: NTuple{3, Int}
    rgrids :: Vector{NPWGrid{Float64}}
    kgrids :: Vector{NPWGrid{Complex{Float64}}}

	FT     :: FFTW.rFFTWPlan{Float64,-1,false,3}
    iFT    :: AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.rFFTWPlan{Complex{Float64},1,false,3},Float64}
end

#==============================================================================#
# Constructors
#==============================================================================#

function FFTBuddy(npw::NTuple{3,<:Integer}; ngrids::Integer = 0, nthreads::Integer = -1)
	tmp_real = zeros(Float64, npw)
	tmp_imag = zeros(Complex{Float64}, floor(Int, npw[1]/2 + 1), npw[2], npw[3])

	# Setup multi-threads and FFT plan
	if nthreads > 0
		FFTW.set_num_threads(nthreads)
	else
		nthreads = determine_nthreads(prod(npw))
		FFTW.set_num_threads(nthreads)
	end
	FT = plan_rfft(tmp_real; flags = FFTW.MEASURE)
	iFT = plan_irfft(tmp_imag, npw[1]; flags = FFTW.MEASURE) #inv(FT)
	fft = FFTBuddy(npw, [], [], FT, iFT)

	# Add the grids
	for i = 1:ngrids
		add_rgrid!(fft)
		add_kgrid!(fft)
	end

	return fft
end

#==============================================================================#
# Methods
#==============================================================================#

show(io::IO, fft::FFTBuddy) = @printf(io, "FFTBuddy(npw = %s)", fft.npw)

num_rgrids(fft::FFTBuddy) = length(fft.rgrids)
num_kgrids(fft::FFTBuddy) = length(fft.kgrids)

add_rgrid!(fft::FFTBuddy) = push!(fft.rgrids, zeros(Float64, fft.npw[1], fft.npw[2], fft.npw[3]))
add_kgrid!(fft::FFTBuddy) = push!(fft.kgrids, zeros(Complex{Float64}, floor(Int, fft.npw[1]/2 + 1), fft.npw[2], fft.npw[3]))

empty_rgrids!(fft::FFTBuddy) = empty!(fft.rgrids)
empty_kgrids!(fft::FFTBuddy) = empty!(fft.kgrids)

function rgrid(fft::FFTBuddy, idx::Integer)
	while num_rgrids(fft) < idx
		add_rgrid!(fft)
	end
	return fft.rgrids[idx]
end

function kgrid(fft::FFTBuddy, idx::Integer)
	while num_kgrids(fft) < idx
		add_kgrid!(fft)
	end
	return fft.kgrids[idx]
end