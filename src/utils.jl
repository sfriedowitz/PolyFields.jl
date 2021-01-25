#==============================================================================#
# Utility functions
#==============================================================================#

ordered_pair(i::Integer, j::Integer) = i <= j ? (i, j) : (j, i)

ndims(dims::NTuple{N,<:Integer}) where {N} = sum(dims .> 1)

ksize(dims::NTuple{3,TI}) where {TI <: Integer} = (floor(TI, dims[1]/2+1), dims[2], dims[3])

function squeeze(A::AbstractArray)
    outshape = Tuple(s for s in size(A) if s > 1)
    return reshape(A, outshape)
end

function determine_nthreads(size::Integer, serial_cut::Integer = 20000)
    return min(Threads.nthreads(), 1 + size÷serial_cut)
end

#==============================================================================#

"""
    copydict!(dest, src)

Copy data from arrays in `src` to `dest`, adding keys to `dest` if necessary.
"""
function copydict!(dest::Dict{TK,TV}, src::Dict{TK,TV}) where {TK,TV<:AbstractArray}
    for (key, val) in src
        if haskey(dest, key) && size(val) == size(dest[key])
            dest[key] .= val
        else
            dest[key] = deepcopy(val)
        end
    end
    return nothing
end

"""
    dotdicts(a, b)

Compute the total dot product between arrays with matching keys in the two dictionaries.
"""
function dotdicts(a::Dict{TK,TV}, b::Dict{TK,TV}) where {TK,TV<:AbstractArray}
    # Performs a dot product between all arrays in a and b dictionaries
    d = 0.0
    for (key, val1) in a
        if haskey(b, key) && size(val1) == size(b[key])
            val2 = b[key]
            d += dot(val1, val2)
        end
    end
    return d
end

"""
    interpolate_grid(grid{TF}, outshape) where {TF <: Real}

Interpolate a real-space field or density grid to the shape defined by `outshape`,
using a cubic b-spline method.

If the interpolation is performed to an output shape of different dimension,
the input grid is either copied or sliced, and the output grid may be inaccurate.
"""
function interpolate_grid(grid::FieldGrid{TF}, outshape::NTuple{3,<:Integer}) where {TF <: Real}
    inshape = size(grid) 

    sqz = squeeze(grid)
    itp = interpolate(sqz, BSpline(Cubic(Line(OnGrid()))))
    itp_dim = ndims(size(itp))

    in_single = Tuple(i for i = 1:3 if inshape[i] == 1)
    out_single = Tuple(i for i = 1:3 if outshape[i] == 1)
    if in_single != out_single
        @warn "Interpolating to different dimensions -- output array may be inaccurate."
    end
    
    # Points at which to evaluate the interpolation function
    hooks = Vector{Float64}[]
    for i = 1:3
        stop_pt = outshape[i] == 1 ? 1 : inshape[i]
        push!(hooks, collect(range(1.0, stop_pt; length = outshape[i])))
    end
    
    # Loop along and interpolate
    itp_grid = zeros(TF, outshape)
    for dx = 1:outshape[1]
        hx = hooks[1][dx]

        for dy = 1:outshape[2]
            hy = hooks[2][dy]

            for dz = 1:outshape[3]
                hz = hooks[3][dz]

                # Choose which dimensions to copy along
                if itp_dim == 3
                    itp_grid[dx, dy, dz] = itp(hx, hy, hz)
                elseif itp_dim == 2
                    itp_grid[dx, dy, dz] = itp(hx, hy)
                else
                    itp_grid[dx, dy, dz] = itp(hx)
                end

            end
        end
    end
    
    return itp_grid
end