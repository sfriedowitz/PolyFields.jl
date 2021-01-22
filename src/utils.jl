#==============================================================================#
# Utility functions
#==============================================================================#

function determine_nthreads(size::Integer, serial_cut::Integer = 20000)
    return min(Threads.nthreads(), 1 + size÷serial_cut)
end

ordered_pair(i::Integer, j::Integer) = i < j ? (i, j) : (j, i)

num_dims(grid::PWGrid) = num_dims(size(grid))

function num_dims(npw::NTuple{3,<:Integer})
    keep_dims = Tuple(i for i in npw if i != 1)
    return length(keep_dims)
end

function squeeze(A::AbstractArray)
     keep_dims = Tuple(i for i in size(A) if i != 1)
     return reshape(A, keep_dims)
end

#==============================================================================#

function copy_dict!(src::Dict{TK,TV}, dst::Dict{TK,TV}) where {TK,TV<:AbstractArray}
    # Copies the contents of dictionary src into dst
    # If key from src does not exist in dst
    #   values of src are deepcopied into dst
    # If key does exist, dst is modified in place
    # Assumes array value sizes are commensurate between src and dst
    for (key, val) in src
        if !haskey(dst, key)
            dst[key] = deepcopy(val)
        else
            dst[key] .= val
        end
    end
    return nothing
end

function dot_dicts(a::Dict{TK,TV}, b::Dict{TK,TV}) where {TK,TV<:AbstractArray}
    # Performs a dot product between all arrays in a and b dictionaries
    d = 0.0
    for (key, val1) in a
        @assert haskey(b, key)
        val2 = b[key]
        d += dot(val1, val2)
    end
    return d
end

#==============================================================================#

"""
    triclinic_vectors(a, b, c, alpha, beta, gamma)

Return a basis matrix for given unit cell side lengths and angles.
"""
function triclinic_vectors(a::Real, b::Real, c::Real, alpha::Real, beta::Real, gamma::Real)
    @assert a > 0 && b > 0 && c > 0
    @assert alpha > 0 && beta > 0 && gamma > 0
    @assert alpha < pi && beta < pi && gamma < pi

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
    triclinic_angles(R_basis)

Return a set of unit cell side lengths and angles for a real basis matrix `R`.
"""
function triclinic_angles(R::Mat3D)
    lx, ly, lz = R[1,1], R[2,2], R[3,3]
    xy, xz, yz = R[1,2], R[1,3], R[2,3]

    a = lx
    b = sqrt(ly^2 + xy^2)
    c = sqrt(lz^2 + xz^2 + yz^2)

    alpha = acos((xy*xz + ly*yz)/(b*c))
    beta  = acos(xz/c)
    gamma = acos(xy/b)

    lengths = @SVector [a, b, c]
    angles = @SVector [alpha, beta, gamma]

    return lengths, angles
end

"""
    interpolate_grid(grid{TF}, outshape) where {TF <: Real}

Interpolate a real-space field or density grid to the shape defined by `outshape`,
using a cubic b-spline method.

If the interpolation is performed to an output shape of different dimension,
the input grid is either copied or sliced, and the output grid may be inaccurate.
"""
function interpolate_grid(grid::PWGrid{TF}, outshape::NTuple{3,<:Integer}) where {TF <: Real}
    inshape = size(grid)            
    grid_sq = squeeze(grid)
    itp = interpolate(grid_sq, BSpline(Cubic(Line(OnGrid()))) )
    itp_dim = length(size(itp))
    
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
    itp_grid = zeros(TF, outshape...)
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