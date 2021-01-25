nmonomers(species::AbstractSpecies) = length(species.mids)

hasmonomer(species::AbstractSpecies, mid::Integer) = mid in species.mids

monomer_fraction(species::AbstractSpecies, mid::Integer) = 0.0

"""
    density!(species)

Update the internal density fields for an `AbstractSpecies`.
"""
density!(species::AbstractSpecies) = nothing

"""
    setup!(species, sys)

Setup internal state of the species for the provided system.
"""
setup!(species::AbstractSpecies, sys::FieldSystem) = nothing

#==============================================================================#

"""
    chain_grid(N_block, Ns)

Calculate the total chain length and block discretization
for a generic polymer chain with variable number of blocks.

The number of contour segments returned for each block needs to be even
for use of Simpson's rule in density/stress calculations.
"""
function chain_grid(N_block::AbstractVector{<:Real}, Ns::Integer)
    N = sum(N_block)
    nblocks = length(N_block)
    ds = N / Ns

    begin_block = Vector{Int}(undef, nblocks+1)
    ds_block = Vector{Float64}(undef, nblocks)

    bgn = 1
    begin_block[1] = bgn

    # Loop over all blocks in the chain
    for iblk = 1:nblocks
        Ns_half = floor(Int, N_block[iblk] / ds / 2.0 + 0.5)
        
        # Calculate step size for each block
        if Ns_half == 0
            Ns_half = 1
            ds_block[iblk] = N_block[iblk] / 2.0
        else
            ds_block[iblk] = N_block[iblk] / Ns_half / 2.0
        end
        
        # Calculate location of block start
        bgn = bgn + Ns_half * 2
        begin_block[iblk+1] = bgn
    end

    return N, begin_block, ds_block
end