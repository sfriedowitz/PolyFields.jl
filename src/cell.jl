"""
    Cell1D(a)
    Cell2D(a, b, gamma = pi/2)
    Cell3D(a, b, c, alpha = pi/2, beta = pi/2, gamma = pi/2)
    Cell3D(lengths, angles = [pi/2, pi/2, pi/2])

A unit cell representing a real-space and inverse-space basis in 1, 2 or 3 dimensions.

The input cell parameters represent:
* `a` : Side length of ``a_1``
* `b` : Side length of ``a_2``
* `c` : Side length of ``a_3``
* `alpha` : Angle between ``a_2`` and ``a_3`` vectors
* `beta`  : Angle between ``a_1`` and ``a_3`` vectors
* `gamma` : Angle between ``a_1`` and ``a_2`` vectors
"""
mutable struct Cell
    dim     :: Int

    # All the fields containing basis information
    params  :: Vector{Float64}
    volume  :: Float64

    R_basis :: SMat3D{Float64} # R[:,i] = Bravais lattice basis vector a_i
    G_basis :: SMat3D{Float64} # G[:,i] = Reciprocal lattice basis vector b_i

    dR      :: Array{Float64,3} # Derivatives of a_i w.r.t. num cell parameters
    dG      :: Array{Float64,3} # Derivatives of b_i
    dRR     :: Array{Float64,3} # Derivatives of a_i ⋅ a_j
    dGG     :: Array{Float64,3} # Derivatives of b_i ⋅ b_j

    # Grids for k-vectors
    ksq     :: PWGrid{Float64}
    dksq    :: Vector{PWGrid{Float64}}

    function Cell(dim::Integer, params::Vector{<:Real})
        cell = new()
        cell.dim = dim
        cell.params = params

        # Setup all of the arrays
        nparam = num_cell_params(cell)
        cell.dR = zeros(3, 3, nparam)
        cell.dG = zeros(3, 3, nparam)
        cell.dRR = zeros(3, 3, nparam)
        cell.dGG = zeros(3, 3, nparam)

        # Add the k-grids
        cell.ksq = PWGrid{Float64}(undef, 0, 0, 0)
        cell.dksq = []

        # Update the basis information to make current before return
        update_cell!(cell)

        return cell
    end
end

#==============================================================================#
# Constructors
#==============================================================================#

function Cell1D(a::Real)
    params = [a]
    return Cell(1, params)
end

function Cell2D(a::Real, b::Real, gamma::Real = pi/2)
    params = [a, b, gamma]
    return Cell(2, params)
end

function Cell3D(a::Real, b::Real, c::Real, alpha::Real = pi/2, beta::Real = pi/2, gamma::Real = pi/2)
    params = [a, b, c, alpha, beta, gamma]
    return Cell(3, params)
end

function Cell3D(lengths::Vector{<:Real}, angles::Vector{<:Real} = [pi/2, pi/2, pi/2])
    @assert length(lengths) == 3 && length(angles) == 3
    return Cell3D(lengths..., angles...)
end

#==============================================================================#
# Methods
#==============================================================================#

show(io::IO, cell::Cell) = @printf(io, "Cell(dim = %d, V = %.2f)", cell.dim, cell.volume)

num_cell_params(cell::Cell) = length(cell.params)

function parameter_string(cell::Cell)
    p = cell.params
    if cell.dim == 1
        return @sprintf "%.3f" p[1]
    elseif cell.dim == 2
        return @sprintf "%.3f, %.3f, %.3f" p[1] p[2] p[3]
    elseif cell.dim == 3
        return @sprintf "%.3f, %.3f, %.3f, %.3f, %.3f, %.3f" p[1] p[2] p[3] p[4] p[5] p[6]
    end
    return ""
end

function setup_kgrids!(cell::Cell, npw::NTuple{3,<:Integer})
    dim = num_dims(npw)
    if dim != cell.dim
        error("Invalid Cell dimension and plane-wave grid dimensions.")
    end

    # Setup the k-squared grid
    cell.ksq = zeros(Float64, floor(Int, npw[1]/2 + 1), npw[2], npw[3])
    ksq_grid!(cell)

    # Setup the grid of k-vectors
    cell.dksq = [zeros(Float64, floor(Int, npw[1]/2 + 1), npw[2], npw[3]) for k = 1:num_cell_params(cell)]
    dksq_grid!(cell)

    return nothing
end

"""
    update_cell!(cell)

Update the cell basis matrices and wavevector grids
given the current set of parameters (lengths/angles).

Must be called after every update to cell parameters during a variable cell calculation.
"""
function update_cell!(cell::Cell)
    # Construct a new basis matrix from cell params
    update_basis!(cell)

    # Calculate cell basis derivatives: dR, dG, dRR, dGG
    d_basis!(cell)
    dd_basis!(cell)

    # Update the ksq and dksq grids
    ksq_grid!(cell)
    dksq_grid!(cell)

    return nothing
end

function update_basis!(cell::Cell)
    p = cell.params
    if cell.dim == 1
        cell.R_basis = triclinic_vectors(p[1], 1.0, 1.0, pi/2, pi/2, pi/2)
    elseif cell.dim == 2
        cell.R_basis = triclinic_vectors(p[1], p[2], 1.0, pi/2, pi/2, p[3])
    else
        cell.R_basis = triclinic_vectors(p...)
    end
    cell.G_basis = 2*pi*transpose(inv(cell.R_basis))
    cell.volume = det(cell.R_basis)
    return nothing
end

function d_basis!(cell::Cell)
    # Calculates the variation of the real and reciprocal basis matrices
    #   w.r.t. each of the (6) cell parameters, a/b/c/alpha/beta/gamma
    # Analytically:
    #   dR_ijk = d(R_ij) / d(theta_k)
    #   dG_ijk = d(G_ij) / d(theta_k)
    fill!(cell.dR, 0.0)
    fill!(cell.dG, 0.0)

    # Calculate the derivatives of the real basis matrix w.r.t each parameter
    #   R == 3x3 matrix
    #   grad(r) = 3x3x6 array, gradient of cell dimensions adds a rank
    if cell.dim == 1
        cell.dR[1,1,1] = 1.0
    elseif cell.dim == 2
        a, b, gamma = cell.params
        cell.dR[1,1,1] = 1.0
        cell.dR[1,2,2] = cos(gamma)
        cell.dR[2,2,2] = sin(gamma)
        cell.dR[1,2,3] = -b * sin(gamma)
        cell.dR[2,2,3] = b * cos(gamma)
    elseif cell.dim == 3
        a, b, c, alpha, beta, gamma = cell.params
        cx, cy, cz = cell.R_basis[1,3], cell.R_basis[2,3], cell.R_basis[3,3]
        cell.dR[1,1,1] = 1.0
        cell.dR[1,2,2] = cos(gamma)
        cell.dR[2,2,2] = sin(gamma)
        cell.dR[1,3,3] = cos(beta)
        cell.dR[2,3,3] = cy / c
        cell.dR[3,3,3] = cz / c
        # The dirty ones w/ angle derivatives, 3rd column gets messy af
        cell.dR[2,3,4] = -c * sin(alpha) / sin(gamma)
        cell.dR[3,3,4] = c^2 * (cos(alpha) - cos(beta)*cos(gamma)) * sin(alpha) / sin(gamma)^2 / cz
        cell.dR[1,3,5] = -c * sin(beta)
        cell.dR[2,3,5] = c * sin(beta) / tan(gamma)
        cell.dR[3,3,5] = c^2 * (-cos(alpha) / tan(gamma) + cos(beta) / sin(gamma)) * sin(beta) / sin(gamma) / cz
        cell.dR[1,2,6] = -b * sin(gamma)
        cell.dR[2,2,6] = b * cos(gamma)
        cell.dR[2,3,6] = c * (-cos(alpha) / tan(gamma) + cos(beta) / sin(gamma)) / sin(gamma)
        cell.dR[3,3,6] = c^2 * ((2 + cos(2*alpha) + cos(2*beta)) / tan(gamma) - cos(alpha) * cos(beta)
             * (3 + cos(2*gamma)) / sin(gamma) )  / sin(gamma)^2 / cz / 2
    end

    # Update the dG matrix using the dR matrix calculated above
    #   d(b_u) = -(1/2pi) * sum_v [ b_u * d(a_v) b_v ] (Jian's thesis & 2003 stress paper)
    #   So we can express reciprocal variation in terms of real variation dR matrix
    G = cell.G_basis
    dR = cell.dR
    for k = 1:num_cell_params(cell)
        for i = 1:cell.dim
            for j = 1:cell.dim
                for l = 1:cell.dim
                    for m = 1:cell.dim
                        cell.dG[i,j,k] -= G[i,l] * dR[m,l,k] * G[m,j]
                    end
                end
            end
        end
    end
    cell.dG ./= 2*pi

    return nothing
end

function dd_basis!(cell::Cell)
    # Calculate derivatives of basis vector dot products w.r.t. cell params
    # Analytically:
    #   dRR_ijk = d(a_i ⋅ a_j) / d(theta_k)
    #   dGG_ijk = d(b_i ⋅ b_j) / d(theta_k)
    fill!(cell.dRR, 0.0)
    fill!(cell.dGG, 0.0)

    # Updates dRR and dGG based on the current basis matrix
    R = transpose(cell.R_basis) # Transpose works to correct indexing from PSCF
    G = cell.G_basis
    dR = cell.dR
    dG = cell.dG
    for k = 1:num_cell_params(cell)
        for i = 1:cell.dim
            for j = 1:cell.dim
                for l = 1:cell.dim
                    cell.dRR[i,j,k] += R[i,l] * dR[l,j,k] + R[j,l] * dR[l,i,k]
                    cell.dGG[i,j,k] += G[l,i] * dG[l,j,k] + G[l,j] * dG[l,i,k]
                end
            end
        end
    end

    return nothing
end

#==============================================================================#

"""
    ksq_grid!(cell)

Update the |k^2| grid given the current reciprocal cell basis `G_basis`.
The wavevectors are shifted in accordance with the output of a real FFT.
"""
function ksq_grid!(cell::Cell)
    Nx, Ny, Nz = size(cell.ksq)
    G = cell.G_basis

    for idx in CartesianIndices(cell.ksq)
        ix, iy, iz = idx[1]-1, idx[2]-1, idx[3]-1

        jx = ix  # Assumes first dim is half-space
        jy = iy <= Ny/2 ? iy : Ny-iy
        jz = iz <= Nz/2 ? iz : Nz-iz 
        jvec = @SVector [jx, jy, jz]
        kvec = G * jvec

        kdot = dot(kvec, kvec)
        @inbounds cell.ksq[idx] = kdot
    end

    return nothing
end

"""
    dksq_grid!(cell)

Update the derivatives of |k^2| for each cell parameter.
"""
function dksq_grid!(cell::Cell)
    Nx, Ny, Nz = size(cell.ksq)
    G = cell.G_basis
    dGG = cell.dGG

    #Threads.@threads for k = 1:length(cell.dksq)
    for k = 1:length(cell.dksq)
        dksq = cell.dksq[k]
        fill!(dksq, 0.0)

        for idx in CartesianIndices(dksq)
            ix, iy, iz = idx[1]-1, idx[2]-1, idx[3]-1

            jx = ix  # Assumes first dim is half-space
            jy = iy <= Ny/2 ? iy : Ny-iy
            jz = iz <= Nz/2 ? iz : Nz-iz 
            jvec = @SVector [jx, jy, jz]

            for l = 1:cell.dim
                for m = 1:cell.dim
                    @inbounds dksq[idx] += jvec[l] * jvec[m] * dGG[l,m,k]
                end
            end

        end
    end

    return nothing
end