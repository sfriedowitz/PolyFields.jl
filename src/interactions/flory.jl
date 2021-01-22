"""
    mutable struct FloryInteraction <: AbstractInteraction

An interaction representing a Flory-Huggins ``\\chi``-interaction
between different monomer types.

The value of the ``\\chi``-parameter is specified with the `set_interaction!` method.
``\\chi``-interactions between like-monomers are not allowed.

### Constructors
```
FloryInteraction()
```
"""
mutable struct FloryInteraction <: AbstractInteraction
    mids   :: Set{Int}
    chis   :: Dict{NTuple{2,Int}, Float64}
    system :: Option{FieldSystem}
    FloryInteraction() = new(Set(), Dict(), nothing)
end

#==============================================================================#
# Methods
#==============================================================================#

show(io::IO, itx::FloryInteraction) = @printf(io, "FloryInteraction(%d monomers, %s pairs)", length(itx.ids), length(itx.chis))

function setup!(itx::FloryInteraction, sys::FieldSystem)
    itx.system = sys
    return nothing
end

function set_interaction!(itx::FloryInteraction, amon::Monomer, bmon::Monomer, chi::Real)
    @assert amon.id != bmon.id

    mids = ordered_pair(amon.id, bmon.id)
    if haskey(itx.chis, mids)
        @warn "Replacing interaction for monomer pair with ids = $(mids)."
    end

    push!(itx.mids, amon.id)
    push!(itx.mids, bmon.id)
    itx.chis[mids] = chi

    return nothing
end

#==============================================================================#

function energy(itx::FloryInteraction)
    @assert !isnothing(itx.system)
    sys = itx.system
    monomers = sys.monomers

    energy = 0.0
    for alpha in itx.mids
        for beta in itx.mids
            pair = ordered_pair(alpha, beta)
            if !haskey(itx.chis, pair); continue; end
            chi = itx.chis[pair]

            va = monomers[alpha].size
            vb = monomers[beta].size
            rho_alpha = sys.density[alpha]
            rho_beta = sys.density[beta]

            @simd for i in eachindex(rho_alpha)
                @inbounds energy += chi * (rho_alpha[i] / va) * (rho_beta[i] / vb)
            end
        end
    end

    return energy / num_grid(sys)
end 

function energy_bulk(itx::FloryInteraction)
    @assert !isnothing(itx.system)
    sys = itx.system
    monomers = sys.monomers

    energy = 0.0
    for alpha in itx.mids
        for beta in itx.mids
            pair = ordered_pair(alpha, beta)
            if !haskey(itx.chis, pair); continue; end
            chi = itx.chis[pair]

            va = monomers[alpha].size
            vb = monomers[beta].size
            rho_alpha = sys.density_bulk[alpha]
            rho_beta = sys.density_bulk[beta]

            energy += chi * (rho_alpha / va) * (rho_beta / vb)
        end
    end

    return energy
end

function add_potential!(itx::FloryInteraction, alpha::Integer, pot::NPWGrid)
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