"""
    mutable struct FloryInteraction <: AbstractInteraction

An interaction representing a Flory-Huggins ``\\chi``-interaction between different monomer types.

The interaction energy has the form:
``H_{int} = (1/2) ∑_i ∑_j χ_{ij} ∫dr ϕ_i(r) ϕ_j(r)``
"""
mutable struct FloryInteraction <: AbstractInteraction
    mids   :: Set{Int}
    chis   :: Dict{NTuple{2,Int}, Float64}
    system :: Option{FieldSystem}
    FloryInteraction() = new(Set(), Dict(), nothing)
end

#==============================================================================#

Base.show(io::IO, itx::FloryInteraction) = @printf(io, "FloryInteraction(%d monomers, %s pairs)", length(itx.mids), length(itx.chis))

function setup!(itx::FloryInteraction, sys::FieldSystem)
    itx.system = sys
    return nothing
end

function set_interaction!(itx::FloryInteraction, alpha::Integer, beta::Integer, chi::Real)
    @assert alpha != beta
    pair = ordered_pair(alpha, beta)
    if haskey(itx.chis, pair)
        @warn "Replacing interaction for monomer pair with IDs = $(pair)."
    end

    push!(itx.mids, alpha)
    push!(itx.mids, beta)
    itx.chis[pair] = chi

    return nothing
end

function energy(itx::FloryInteraction)
    @assert !isnothing(itx.system)
    sys = itx.system

    energy = 0.0
    for (pair, chi) in itx.chis
        alpha, beta = pair
        if hasmonomer(sys, alpha) && hasmonomer(sys, beta)
            va = sys.monomers[alpha].vol
            vb = sys.monomers[beta].vol
            rho_alpha = sys.density[alpha]
            rho_beta = sys.density[beta]

            @simd for i in eachindex(rho_alpha)
                @inbounds energy += chi * (rho_alpha[i] / va) * (rho_beta[i] / vb)
            end
        end
    end

    return energy / ngrid(sys)
end 

function energy_bulk(itx::FloryInteraction)
    @assert !isnothing(itx.system)
    sys = itx.system
    monomers = sys.monomers

    energy = 0.0
    for (pair, chi) in itx.chis
        alpha, beta = pair
        if hasmonomer(sys, alpha) && hasmonomer(sys, beta)
            va = monomers[alpha].vol
            vb = monomers[beta].vol
            rho_alpha = sys.monomer_fracs[alpha]
            rho_beta = sys.monomer_fracs[beta]

            energy += chi * (rho_alpha / va) * (rho_beta / vb)
        end
    end

    return energy
end

function potential!(itx::FloryInteraction, alpha::Integer, pot::FieldGrid)
    @assert !isnothing(itx.system)
    sys = itx.system

    for beta in itx.mids
        pair = ordered_pair(alpha, beta)
        if !haskey(itx.chis, pair); continue; end
        chi = itx.chis[pair]

        rho_beta = sys.density[beta]
        vb = sys.monomers[beta].vol

        @simd for i in eachindex(pot)
            @inbounds pot[i] += chi * (rho_beta[i] / vb)
        end
    end

    return nothing
end