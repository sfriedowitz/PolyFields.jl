#==============================================================================#
# Utility functions
#==============================================================================#

ndims(dims::NTuple{N,<:Integer}) where {N} = sum(dims .> 1)

ksize(dims::NTuple{3,TI}) where {TI <: Integer} = (floor(TI, dims[1]/2+1), dims[2], dims[3])

function squeeze(A::AbstractArray)
    outshape = Tuple(s for s in size(A) if s > 1)
    return reshape(A, outshape)
end

ordered_pair(i::Integer, j::Integer) = i <= j ? (i, j) : (j, i)

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

#==============================================================================#

"""
    make_basis(a, b, c, alpha, beta, gamma)

Return a basis matrix for given unit cell side lengths and angles.
"""
function make_basis(a::Real, b::Real, c::Real, alpha::Real, beta::Real, gamma::Real)
    p2 = pi/2
    if isapprox(alpha,p2) && isapprox(beta,p2) && isapprox(gamma,p2)
        # We have an orthogonal box
        return SMat3D(a, 0.0, 0.0, 0.0, b, 0.0, 0.0, 0.0, c)
    else
        # We have a triclinic box
        # Use exact trigonometric values for right angles
        if isapprox(alpha, p2)
            cos_alpha = 0.0
        else
            cos_alpha = cos(alpha)
        end
        if isapprox(beta, p2)
            cos_beta = 0.0
        else
            cos_beta = cos(beta)
        end
        if isapprox(gamma, p2)
            cos_gamma = 0.0
            sin_gamma = 1.0
        else
            cos_gamma = cos(gamma)
            sin_gamma = sin(gamma)
        end

        # Convert from angle representation to basis components
        ax = a
        bx = b * cos_gamma
        by = b * sin_gamma
        cx = c * cos_beta
        cy = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
        cz = sqrt(c^2 - cx^2 - cy^2)

        return SMat3D(ax, 0.0, 0.0, bx, by, 0.0, cx, cy, cz)
    end
end

"""
    basis_dimensions(basis)

Return a vector of unit cell side lengths and angles for a basis matrix in real space.
"""
function basis_dimensions(basis::Mat3D)
    lx, ly, lz = basis[1,1], basis[2,2], basis[3,3]
    xy, xz, yz = basis[1,2], basis[1,3], basis[2,3]

    a = lx
    b = sqrt(ly^2 + xy^2)
    c = sqrt(lz^2 + xz^2 + yz^2)

    alpha = acos((xy*xz + ly*yz)/(b*c))
    beta  = acos(xz/c)
    gamma = acos(xy/b)

    return @SVector [a, b, c, alpha, beta, gamma]
end

#==============================================================================#

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