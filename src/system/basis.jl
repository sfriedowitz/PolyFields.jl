# Recognized crystal types

const CRYSTAL1D = (:lamellar,)
const CRYSTAL2D = (:square, :rectangular, :hexagonal, :oblique)
const CRYSTAL3D = (:cubic, :tetragonal, :orthorhombic, :hexagonal, :monoclinic, :triclinic, :pnna, :fddd1)

#==============================================================================#

make_basis(dim::Integer, crystal::Symbol, params...) = make_basis(dim, crystal, [params...])

function make_basis(dim::Integer, crystal::Symbol, params::AbstractVector{<:Real})
    basis = @SMatrix zeros(3,3)

    if dim == 1
        if !(crystal in CRYSTAL1D)
            error("Unrecognized crystal system `$(crystal)` for dim = $(dim).")
        end
        if crystal == :lamellar
            @assert length(params) == 1
            basis = basis_matrix(params[1], 1.0, 1.0, pi/2, pi/2, pi/2)
        end
    elseif dim == 2
        if !(crystal in CRYSTAL2D)
            error("Unrecognized crystal system `$(crystal)` for dim = $(dim).")
        end
        if crystal == :square
            @assert length(params) == 1
            basis = basis_matrix(params[1], params[1], 1.0, pi/2, pi/2, pi/2)
        elseif crystal == :rectangular
            @assert length(params) == 2
            basis = basis_matrix(params[1], params[2], 1.0, pi/2, pi/2, pi/2)
        elseif crystal == :hexagonal
            @assert length(params) == 1
            basis = basis_matrix(params[1], params[1], 1.0, pi/2, pi/2, pi/3)
        elseif crystal == :oblique
            @assert length(params) == 3
            basis = basis_matrix(params[1], params[2], 1.0, pi/2, pi/2, params[3])
        end
    elseif dim == 3
        if !(crystal in CRYSTAL3D)
            error("Unrecognized crystal system `$(crystal)` for dim = $(dim).")
        end
        if crystal == :cubic
            @assert length(params) == 1
            basis = basis_matrix(params[1], params[1], params[1], pi/2, pi/2, pi/2)
        elseif crystal == :tetragonal
            @assert length(params) == 2
            basis = basis_matrix(params[1], params[1], params[2], pi/2, pi/2, pi/2)
        elseif crystal == :orthorhombic
            @assert length(params) == 3
            basis = basis_matrix(params[1], params[2], params[3], pi/2, pi/2, pi/2)
        elseif crystal == :hexagonal
            @assert length(params) == 2
            basis = basis_matrix(params[1], params[1], params[2], pi/2, pi/2, pi/3)
        elseif crystal == :monoclinic
            @assert length(params) == 4
            basis = basis_matrix(params[1], params[2], params[3], pi/2, params[4], pi/2)
        elseif crystal == :triclinic
            @assert length(params) == 6
            basis = basis_matrix(params...)
        elseif crystal == :pnna
            @assert length(params) == 1
            basis = basis_matrix(2.0*params[1], sqrt(3.0)*params[1], params[1], pi/2, pi/2, pi/2)
        elseif crystal == :fddd1
            @assert length(params) == 1
            basis = basis_matrix(params[1], sqrt(3.0)*params[1], 2.0*sqrt(3.0)*params[1], pi/2, pi/2, pi/2)
        end
    end

    return basis
end

#==============================================================================#

"""
    basis_matrix(a, b, c, alpha, beta, gamma)

Return a basis matrix for given unit cell side lengths and angles.
"""
function basis_matrix(a::Real, b::Real, c::Real, alpha::Real = pi/2, beta::Real = pi/2, gamma::Real = pi/2)
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