#==============================================================================#
# Field IO
#==============================================================================#

"""
"""
function save_fields(sys::FieldSystem, fpath::AbstractString; type::AbstractString = "omega")    
    # Open file handle
    io = open(fpath, "w")
    
    # For ordering output columns
    mids = sort(collect(keys(sys.fields)))
    
    # Headers
    version = "1"
    field_type = "R" # Hard-coded for now
    dims = dimensions(sys.cell)
    npw = sys.npw
    
    @printf(io, "# Generated by PolyFields.jl on %s\n", Dates.now()) 
    @printf(io, "# Version: %s\n", version)
    print(io, "#\n")
    
    print(io,   "# System:\n")
    @printf(io, "# NPW  = %s\n", npw)
    @printf(io, "# CELL = [%.3f, %.3f, %.3f, %.3f, %.3f, %.3f]\n", dims...)
    @printf(io, "# IDS  = %s\n", mids)
    @printf(io, "# TYPE = %s\n", field_type)
    print(io, "#\n")
    
    # Loop over and write the field output
    @printf(io, "# Columns: ix iy iz fields\n")
    for iz = 1:npw[3]
        for iy = 1:npw[2]
            for ix = 1:npw[1]
                @printf(io, "%5d %5d %5d", ix, iy, iz)
                for mid in mids
                    omega = sys.fields[mid]
                    @printf(io, "%15e", omega[ix, iy, iz])
                end
                print(io, "\n")
            end
            if npw[1] > 1
                print(io, "\n")
            end
        end
        if npw[2] > 1
            print(io, "\n")
        end
    end

    close(io)
    return nothing
end

"""
"""
function load_fields(fpath::AbstractString)
    # Right now the line numbers are hard coded 
    # Can we write a more general input function that searches for appropriate code blocks?
    # Should we switch to a structured data format to allow easier reading/searching?

    farr = readlines(fpath)
    if parse(Int, farr[2][end]) != 1
        error("Unknown field format version.")
    end

    # Parse system information
    sys_lines = farr[5:8]

    npw = Tuple(parse(Int, N) for N in split(sys_lines[1][11:end-1], ","))
    dims = parse.(Float64, split(sys_lines[2][11:end-1], ","))
    mids = parse.(Int, split(sys_lines[3][11:end-1], ","))

    field_type = string(sys_lines[4][10])

    # Create dictionary of fields
    fields = Dict{Int, PWGrid{Float64}}()
    for mid in mids
        if field_type == "R"
            fields[mid] = zeros(Float64, npw)
        elseif field_type == "K"
            fields[mid] = zeros(Complex{Float64}, npw)
        end
    end

    # Loop over lines and load field entries
    for line in farr[11:end]
        if !isempty(line)
            vals = split(line)
            ix = parse(Int, vals[1])
            iy = parse(Int, vals[2])
            iz = parse(Int, vals[3])

            # Loop over the rest of the fields
            for (midx, col) in enumerate(4:length(vals))
                v = vals[col]
                mid = mids[midx]
                fields[mid][ix, iy, iz] = parse(Float64, v)
            end
        end
    end

    return fields
end