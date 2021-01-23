"""
    mutable struct GaussianEdwards <: AbstractInteraction

An analytical form of Gaussian fluctuation electrostatics 
using Edward's approximation for all polymer form factors.

### Constructors
```
GaussianEdwards(; lB = 2.287)
```
"""
mutable struct GaussianEdwards <: AbstractInteraction
    lB     :: Float64
    pref   :: Float64
    cphi   :: Float64

    pols   :: Vector{Multiblock}
    salts  :: Vector{Point}

    PP     :: FieldGrid{Float64}
    PS     :: FieldGrid{Float64}
    DPP    :: FieldGrid{Float64}
    DPS    :: FieldGrid{Float64}
    system :: Option{FieldSystem}

    function GaussianEdwards(; lB::Real = 2.287)
        pref = (4*pi*lB)^(3/2) / (12*pi)
        cphi = sqrt(3/(lB*pi))
        return new(lB, pref, cphi, [], [], [], [], [], [], nothing)
    end
end

#==============================================================================#
# Methods
#==============================================================================#

Base.show(io::IO, itx::GaussianEdwards) = @printf(io, "GaussianEdwards(lB = %.3g)", itx.lB)

function setup!(itx::GaussianEdwards, sys::FieldSystem)
    itx.system = sys
    itx.PP = zeros(sys.npw)
    itx.PS = zeros(sys.npw)
    itx.DPP = zeros(sys.npw)
    itx.DPS = zeros(sys.npw)
    return nothing
end

function set_interaction!(itx::GaussianEdwards, species::Multiblock)
    push!(itx.pols, species)
    return nothing
end

function set_interaction!(itx::GaussianEdwards, species::Point)
    push!(itx.salts, species)
    return nothing
end

#==============================================================================#

function energy(itx::GaussianEdwards)
    @assert !isnothing(itx.system)
    sys = itx.system

    # Prefactors and constants
    pref = itx.pref
    cphi = itx.cphi

    energy = 0.0
    for i = 1:num_grid(sys) # Our summation is over all discretized points of the density fields
        # Polymer contributions
        # PP = sum [ (sig^2 phiP / omega / b^2)^(1/2) ]
        PP = 0.0
        for pol in itx.pols
            for (iblk, mid) in enumerate(pol.mid_block)
                bblk = pol.b_block[iblk]
                mon = sys.monomers[mid]
                rho = sys.density[mid]

                PP += sqrt(mon.charge^2 * rho[i] / bblk^2 / mon.size)
            end

        end

        # Salt contributions
        # PS = sum [ (phiS / omega) ]
        PS = 0.0
        for salt in itx.salts
            mid = salt.mids[1]
            mon = sys.monomers[mid]
            rho = sys.density[mid]

            PS += abs(mon.charge) * rho[i] / mon.size # Do we add the salt monomer charge here for multivalent?
        end

        # Free energy density evaluated at local grid point
        fcorr = pref * (cphi * PP - PS) * (2 * cphi * PP + PS)^(1/2)
        energy += fcorr
    end

    return energy / num_grid(sys)
end 

function energy_bulk(itx::GaussianEdwards)
    @assert !isnothing(itx.system)
    sys = itx.system
    monomers = sys.monomers

    energy = 0.0

    return energy
end

function add_potential!(itx::GaussianEdwards, alpha::Integer, pot::FieldGrid)
    @assert !isnothing(itx.system)
    sys = itx.system
    monomers = sys.monomers

    for beta in itx.mids
        pair = ordered_pair(alpha, beta)
        if !haskey(itx.chis, pair); continue; end
        chi = itx.chis[pair]

        vb = monomers[beta].size
        rho_beta = sys.density[beta]

        @simd for i in eachindex(pot)
            @inbounds pot[i] += chi * (rho_beta[i] / vb)
        end
    end

    return nothing
end

















