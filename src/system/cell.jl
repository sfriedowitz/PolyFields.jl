"""
    mutable struct UnitCell

    Cell(dim, crystal, params)

A unit cell representing a real-space and inverse-space basis in 1, 2 or 3 dimensions.

A specified crystal system controls the geometry of the cell. These include
*1D: (:lamellar)
*2D: (:square, :rectangular, :hexagonal, :oblique)
*3D: 
"""
mutable struct UnitCell
    dim     :: Int
    crystal :: Symbol
    volume  :: Float64
    params  :: Vector{Float64}

    R       :: SMat3D{Float64}  # R[:,i] = Bravais lattice basis vector a_i
    G       :: SMat3D{Float64}  # G[:,i] = Reciprocal lattice basis vector b_i
    dR      :: Array{Float64,3} # Derivatives of a_i w.r.t. cell parameters
    dG      :: Array{Float64,3} # Derivatives of b_i
    dRR     :: Array{Float64,3} # Derivatives of a_i ⋅ a_j
    dGG     :: Array{Float64,3} # Derivatives of b_i ⋅ b_j

    # Grids for k-vectors
    ksq     :: FieldGrid{Float64}
    dksq    :: Vector{FieldGrid{Float64}}

    function UnitCell(dim::Integer, crystal::Symbol, params::AbstractVector{<:Real})
        cell = new()
        cell.dim = dim
        cell.crystal = crystal
        cell.params = params

        npr = nparams(cell)
        cell.dR = zeros(3, 3, npr)
        cell.dG = zeros(3, 3, npr)
        cell.dRR = zeros(3, 3, npr)
        cell.dGG = zeros(3, 3, npr)

        # Add the k-grids
        cell.ksq = FieldGrid{Float64}(undef, 0, 0, 0)
        cell.dksq = []

        # Update the basis information to make current before return
        update!(cell)

        return cell
    end
end

UnitCell(dim::Integer, crystal::Symbol, params...) = UnitCell(dim, crystal, [params...])

#==============================================================================#

Base.show(io::IO, cell::UnitCell) = @printf(io, "Cell(dim = %d, V = %.2f)", cell.dim, cell.volume)

nparams(cell::UnitCell) = length(cell.params)

"""
    setup!(cell, sys)

Allocate k-grids of an appropriate size for the system.
"""
function setup!(cell::UnitCell, sys::AbstractSystem)
    dim = ndims(sys.dims)
    if dim != cell.dim
        error("Number of non-singleton dimensions in System and Cell do not match!")
    end

    # Setup the cell
    cell.ksq = zeros(Float64, ksize(sys.dims))
    cell.dksq = [zeros(Float64, ksize(sys.dims)) for k = 1:nparams(cell)]
    
    ksq!(cell)
    dksq!(cell)

    return nothing
end

"""
    update!(cell)

Update the cell basis matrices and wavevector grids
given the current set of parameters (lengths/angles).

Must be called after every update to cell parameters during a variable cell calculation.
"""
function update!(cell::UnitCell)
    # Construct a new basis matrix from cell params
    basis!(cell)

    # Calculate cell basis derivatives: dR, dG, dRR, dGG
    dbasis!(cell)
    ddbasis!(cell)

    # Update the ksq and dksq grids
    ksq!(cell)
    dksq!(cell)

    return nothing
end

#==============================================================================#

function basis!(cell::UnitCell)
    cell.R = make_basis(cell.dim, cell.crystal, cell.params)
    cell.G = 2*pi*transpose(inv(cell.R))
    cell.volume = det(cell.R)
    return nothing
end

function dbasis!(cell::UnitCell)
    # Calculates the variation of the real and reciprocal basis matrices
    #   w.r.t. each of the (6) cell parameters, a/b/c/alpha/beta/gamma
    # Analytically:
    #   dR_ijk = d(R_ij) / d(theta_k)
    #   dG_ijk = d(G_ij) / d(theta_k)
    fill!(cell.dR, 0.0)
    fill!(cell.dG, 0.0)

    # Calculate the derivatives of the real basis matrix w.r.t each parameter
    #   R == 3x3 matrix
    #   grad(R) = 3x3x6 array, gradient of cell dimensions adds a rank
    if cell.dim == 1
        if cell.crystal == :lamellar
            cell.dR[1,1,1] = 1.0
        end
    elseif cell.dim == 2
        if cell.crystal == :square
            cell.dR[1,1,1] = 1.0
            cell.dR[2,2,1] = 1.0
        elseif cell.crystal == :rectangular
            cell.dR[1,1,1] = 1.0
            cell.dR[2,2,2] = 1.0
        elseif cell.crystal == :hexagonal
            cell.dR[1,1,1] = 1.0
            cell.dR[1,2,1] = 0.5
            cell.dR[2,2,1] = sqrt(0.75)
        elseif cell.crystal == :oblique
            a, b, gamma = cell.params
            cell.dR[1,1,1] = 1.0
            cell.dR[1,2,2] = cos(gamma)
            cell.dR[2,2,2] = sin(gamma)
            cell.dR[1,2,3] = -b * sin(gamma)
            cell.dR[2,2,3] = b * cos(gamma)
        end
    elseif cell.dim == 3
        if cell.crystal == :cubic
            cell.dR[1,1,1] = 1.0
            cell.dR[2,2,1] = 1.0
            cell.dR[3,3,1] = 1.0
        elseif cell.crystal == :tetragonal
            cell.dR[1,1,1] = 1.0
            cell.dR[2,2,1] = 1.0
            cell.dR[3,3,2] = 1.0
        elseif cell.crystal == :orthorhombic
            cell.dR[1,1,1] = 1.0
            cell.dR[2,2,2] = 1.0
            cell.dR[3,3,3] = 1.0
        elseif cell.crystal == :hexagonal
            cell.dR[1,1,1] = 1.0
            cell.dR[1,2,1] = 0.5
            cell.dR[2,2,1] = sqrt(0.75)
            cell.dR[3,3,2] = 1.0
        elseif cell.crystal == :monoclinic
            a, b, c, beta = cell.params
            cell.dR[1,1,1] = 1.0
            cell.dR[2,2,2] = 1.0
            cell.dR[1,3,3] = cos(beta)
            cell.dR[3,3,3] = sin(beta)
            cell.dR[1,3,4] = -c*sin(beta)
            cell.dR[3,3,4] = c*cos(beta)
        elseif cell.crystal == :triclinic
            a, b, c, alpha, beta, gamma = cell.params
            cx, cy, cz = cell.R[1,3], cell.R[2,3], cell.R[3,3]
            cell.dR[1,1,1] = 1.0
            cell.dR[1,2,2] = cos(gamma)
            cell.dR[2,2,2] = sin(gamma)
            cell.dR[1,3,3] = cos(beta)
            cell.dR[2,3,3] = cy / c
            cell.dR[3,3,3] = cz / c
            # The dirty ones w/ angle derivatives, 3rd column gets messy
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
        elseif cell.crystal == :pnna
            cell.dR[1,1,1] = 2.0
            cell.dR[2,2,1] = sqrt(3.0)
            cell.dR[3,3,1] = 1.0
        elseif cell.crystal == :fddd1
            cell.dR[1,1,1] = 1.0
            cell.dR[2,2,1] = sqrt(3.0)
            cell.dR[3,3,1] = 2.0*sqrt(3.0)
        end
    end

    # Update the dG matrix using the dR matrix calculated above
    #   d(b_u) = -(1/2pi) * sum_v [ b_u * d(a_v) b_v ] (Jian's thesis & 2003 stress paper)
    #   So we can express reciprocal variation in terms of real variation dR matrix
    G = cell.G
    dR = cell.dR
    for k = 1:nparams(cell)
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

function ddbasis!(cell::UnitCell)
    # Calculate derivatives of basis vector dot products w.r.t. cell params
    # Analytically:
    #   dRR_ijk = d(a_i ⋅ a_j) / d(theta_k)
    #   dGG_ijk = d(b_i ⋅ b_j) / d(theta_k)
    fill!(cell.dRR, 0.0)
    fill!(cell.dGG, 0.0)

    # Updates dRR and dGG based on the current basis matrix
    R = transpose(cell.R) # Transpose works to correct indexing from PSCF
    G = cell.G
    dR = cell.dR
    dG = cell.dG
    for k = 1:nparams(cell)
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
    ksq!(cell)

Update the |k^2| grid given the current reciprocal cell basis `Gbasis`.
The wavevectors are shifted in accordance with the output of a real FFT.
"""
function ksq!(cell::UnitCell)
    Nx, Ny, Nz = size(cell.ksq)
    G = cell.G

    for idx in CartesianIndices(cell.ksq)
        ix, iy, iz = idx[1]-1, idx[2]-1, idx[3]-1

        jx = ix  # Assumes first dim is half-space
        jy = iy <= Ny/2 ? iy : Ny-iy
        jz = iz <= Nz/2 ? iz : Nz-iz 
        jvec = @SVector [jx, jy, jz]
        kvec = G * jvec
        
        @inbounds cell.ksq[idx] = dot(kvec, kvec)
    end

    return nothing
end

"""
    dksq!(cell)

Update the derivatives of |k^2| for each cell parameter.
"""
function dksq!(cell::UnitCell)
    Nx, Ny, Nz = size(cell.ksq)
    G = cell.G
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