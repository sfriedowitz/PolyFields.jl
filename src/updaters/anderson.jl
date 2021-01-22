"""
    AndersonUpdater(; nhist = 50, iter_euler = 1, lam_euler = 0.01)

A field updating scheme implementing Anderson mixing, which obtains a new field configuration
based on a combination of previous configurations and residuals.

The number of configurations stored in history is set by the `nhist` argument.
A given number `iter_euler` of simple Euler update steps with size `lam_euler` 
can be taken at the beginning of iteration to build of the necessary configuration histories.
"""
mutable struct AndersonUpdater <: AbstractFieldUpdater
    nhist      :: Int
    iter       :: Int
    iter_hist  :: Int
    iter_euler :: Int
    lam        :: Float64
    lam_euler  :: Float64
    setup      :: Bool

    res_hist   :: Vector{Dict{Int,NPWGrid{Float64}}} # Residual history for each field
    field_hist :: Vector{Dict{Int,NPWGrid{Float64}}} # Omega field history

    # Coefficient vectors, allocated at each step
    Umn        :: Matrix{Float64}
    Vm         :: Vector{Float64}
    Cn         :: Vector{Float64}

    # Step directions
    WW         :: NPWGrid{Float64}
    DD         :: NPWGrid{Float64}

    function AndersonUpdater(; nhist::Integer = 50, iter_euler::Integer = 1, lam_euler::Real = 0.01)
        # Must take at least 1 initial step before mixing
        if iter_euler <= 1; iter_euler = 1; end

        res_hist = [Dict() for i = 1:nhist+1]
        field_hist = [Dict() for i = 1:nhist+1]

        return new(nhist, 0, 0, iter_euler, 0.0, lam_euler, false, res_hist, field_hist, 
            zeros(nhist, nhist), zeros(nhist), zeros(nhist),
            zeros(0, 0, 0), zeros(0, 0, 0)
        )
    end
end

#==============================================================================#
# Methods
#==============================================================================#

show(io::IO, updater::AndersonUpdater) = @printf(io, "AndersonUpdater(nhist = %d, iter = %d)", updater.nhist, updater.iter)

function setup!(updater::AndersonUpdater, sys::FieldSystem)
    # If already bound to same system and setup, keep histories
    if updater.setup
        @warn "Using previous iterations in AndersonUpdater history."
    else
        if sys.ensemble != Canonical
            error("Anderson mixing can only be used for systems in the canonical ensemble.")
        end

        # Allocate step directions
        updater.WW = zeros(sys.npw)
        updater.DD = zeros(sys.npw)

        # Allocate history arrays
        for dict in updater.res_hist
            for (mid, res) in sys.residuals
                dict[mid] = zeros(sys.npw)
            end
        end
        for dict in updater.field_hist
            for (mid, omega) in sys.fields
                dict[mid] = zeros(sys.npw)
            end
        end

        # Store initial fields
        for (mid, omega) in sys.fields
            updater.field_hist[1][mid] .= omega
        end

        # Reset fields
        updater.iter = 0
        updater.iter_hist = 0
        updater.lam = 0.0
        updater.setup = true
    end

    return nothing
end

#==============================================================================#

function update_state!(sys::FieldSystem, updater::AndersonUpdater)
    # Calculate lam and iter_hist
    updater.iter += 1
    if updater.iter < updater.nhist+1
        updater.iter_hist = updater.iter - 1
        updater.lam = 1.0 - 0.9^updater.iter
    else
        updater.iter_hist = updater.nhist
        updater.lam = 1.0
    end

    # Iteration numbers
    iter = updater.iter
    iter_hist = updater.iter_hist

    # Update density and residuals
    make_residuals!(sys)

    # Update the residual history
    for (mid, res) in sys.residuals
        updater.res_hist[iter_hist+1][mid] .= res
    end

    # First history iteration
    if iter <= updater.iter_euler
        # Update with simple Euler step
        for (mid, omega) in sys.fields
            @. omega += updater.lam_euler * sys.residuals[mid]
        end
    # Second iteration and beyond
    else
        # Reallocate appropriately sized coefficient arrays if needed
        if iter < updater.nhist+2
            updater.Umn = zeros(iter_hist, iter_hist)
            updater.Vm = zeros(iter_hist)
            updater.Cn = zeros(iter_hist)
        else
            fill!(updater.Umn, 0.0)
            fill!(updater.Vm, 0.0)
            fill!(updater.Cn, 0.0)
        end

        # Calculate coefficients of Umn and Vm arrays
        res_iter = updater.res_hist[iter_hist+1]
        field_iter = updater.field_hist[iter_hist+1]

        for mm = 1:iter_hist
            res_mm = updater.res_hist[iter_hist+1-mm]

            for nn = mm:iter_hist
                res_nn = updater.res_hist[iter_hist+1-nn]

                # Compute difference coefficients of Umn matrix
                #   elem = (curr - mm) * (curr - nn)
                for (mid, res) in res_iter
                    rmm = res_mm[mid]
                    rnn = res_nn[mid]

                    elem = 0.0
                    @simd for i in eachindex(res)
                        @inbounds elem += (res[i] - rmm[i]) * (res[i] - rnn[i])
                    end

                    updater.Umn[mm, nn] += elem
                end
                updater.Umn[nn, mm] = updater.Umn[mm, nn]
            end

            # Compute coefficients of Vm matrix
            #   elem = (curr - mm) * curr
            for (mid, res) in res_iter
                rmm = res_mm[mid]

                elem = 0.0
                @simd for i in eachindex(res)
                    @inbounds elem += (res[i] - rmm[i]) * res[i]
                end

                updater.Vm[mm] += elem
            end
        end

        # Inverting Umn matrix
        try
            inv!(lu!(updater.Umn))
            mul!(updater.Cn, updater.Umn, updater.Vm)
        catch err
            error("Error inverting Umn matrix on iteration $(iter): $(err)")
        end

        # Update the fields
        for (mid, omega) in sys.fields
            updater.WW .= omega
            updater.DD .= res_iter[mid]

            for nn = 1:iter_hist
                res_nn = updater.res_hist[iter_hist+1-nn]
                field_nn = updater.field_hist[iter_hist+1-nn]

                @. updater.WW += updater.Cn[nn] * (field_nn[mid] - field_iter[mid])
                @. updater.DD += updater.Cn[nn] * (res_nn[mid] - res_iter[mid])
            end

            @. omega = updater.WW + updater.lam * updater.DD
        end        
    end

    # Update the history arrays
    if iter < updater.nhist+1
        # If not full, add to the array end
        for (mid, omega) in sys.fields
            updater.field_hist[iter+1][mid] .= omega
        end
    else
        # If full, shift all down one slot
        for nh = 1:updater.nhist
            for mid in keys(sys.fields)
                updater.res_hist[nh][mid] .= updater.res_hist[nh+1][mid]
                updater.field_hist[nh][mid] .= updater.field_hist[nh+1][mid]
            end
        end
        for (mid, omega) in sys.fields
            updater.field_hist[end][mid] .= omega
        end
    end

    return nothing
end