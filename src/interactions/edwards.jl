"""
    mutable struct EdwardsInteraction <: AbstractInteraction

An interaction representing an Edward's excluded volume interaction between like-species.

The interaction energy is of the form:
``H = (1/2) ∑_i u_i ∫dr ϕ_i^2(r)``
"""
mutable struct EdwardsInteraction <: AbstractInteraction
    mids   :: Set{Int}
    uvals  :: Dict{Int,Float64}
    system :: Option{FieldSystem}
    EdwardsInteraction() = new(Set(), Dict(), nothing)
end

#==============================================================================#

Base.show(io::IO, itx::EdwardsInteraction) = @printf(io, "EdwardsInteraction(%d monomers)", length(itx.mids))

function setup!(itx::EdwardsInteraction, sys::FieldSystem)
    itx.system = sys
    return nothing
end

function set_interaction!(itx::EdwardsInteraction, mid::Integer, u::Real)
    if mid in itx.mids
        @warn "Replacing interaction for monomer with id = $(mid)."
    end
    push!(itx.mids, mid)
    itx.uvals[mid] = u
    return nothing
end

function energy(itx::EdwardsInteraction)
    @assert !isnothing(itx.system)
    sys = itx.system

    energy = 0.0
    for mid in itx.mids
        if hasmonomer(sys, mid)
            u = itx.uvals[mid]
            v = sys.monomers[mid].vol
            rho = sys.density[mid]

            @simd for i in eachindex(rho)
                @inbounds energy += (u/2)*(rho[i]/v)^2
            end
        end
    end

    return energy/ngrid(sys)
end 

function energy_bulk(itx::EdwardsInteraction)
    @assert !isnothing(itx.system)
    sys = itx.system

    energy = 0.0
    for mid in itx.mids
        if hasmonomer(sys, mid)
            u = itx.uvals[mid]
            v = sys.monomers[mid].vol
            rho = sys.monomer_fracs[mid]

            energy += (u/2)*(rho/v)^2
        end
    end

    return energy
end

function potential!(itx::EdwardsInteraction, mid::Integer, pot::FieldGrid)
    @assert !isnothing(itx.system)
    sys = itx.system

    if mid in itx.mids
        u = itx.uvals[mid]
        v = sys.monomers[mid].vol
        rho = sys.density[mid]

        @simd for i in eachindex(pot)
            @inbounds pot[i] += u*(rho[i]/v)
        end
    end

    return nothing
end