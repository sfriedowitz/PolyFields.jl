"""
    mutable struct EdwardsInteraction <: AbstractInteraction

An interaction representing an Edward's excluded volume interaction between like-species.

The interaction energy is of the form:
``H = ∑_i u_i ∫dr ϕ_i^2(r)``
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

function set_interaction(itx::EdwardsInteraction, mid::Integer, u::Real)
    if mid in itx.mids
        @warn "Replacing interaction for monomer with id = $(mid)."
    end
    itx.uvals[mid] = u
    return nothing
end

function energy(itx::EdwardsInteraction)
    @assert !isnothing(itx.system)
    sys = itx.system
    monomers = sys.monomers

    energy = 0.0
    for alpha in itx.mids
        u = itx.uvals[alpha]
        v = monomers[alpha].size

        rho = sys.density[alpha]
        @simd for i in eachindex(rho)
            @inbounds energy += (u/2)*(rho[i]/v)^2
        end
    end

    return energy/ngrid(sys)
end 

function energy_bulk(itx::EdwardsInteraction)
    @assert !isnothing(itx.system)
    sys = itx.system
    monomers = sys.monomers

    energy = 0.0
    for mid in itx.mids
        u = itx.uvals[mid]
        v = monomers[mid].size

        rho = sys.density_bulk[mid]
        energy += (u/2)*(rho/v)^2
    end

    return energy
end

function potential!(itx::EdwardsInteraction, mid::Integer, pot::FieldGrid)
    @assert !isnothing(itx.system)
    sys = itx.system
    monomers = sys.monomers

    rho = sys.density[mid]
    u = itx.uvals[mid]
    v = monomers[mid].size

    @simd for i in eachindex(pot)
        @inbounds pot[i] += u*(rho[i]/v)
    end

    return nothing
end