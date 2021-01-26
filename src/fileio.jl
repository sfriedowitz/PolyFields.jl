#==============================================================================#
# Field initialization
#==============================================================================#

"""
    fieldinit!(sys; seed = -1, scale = 1.0)
    fieldinit!(sys, fields; field_type = :omega, interpolate = true)
    fieldinit!(sys, fpath; load_cell = :true, interpolate = true)

Initialize fields for all monomer types in the system.
Fields can be initialized randomly with no provided data,
loaded from a field dictionary, or loaded from file.
"""
function fieldinit!(sys::FieldSystem; seed::Integer = -1, scale::Real = 0.1)
    if seed >= 0
        Random.seed!(seed)
    end
    for (mid, omega) in sys.fields
        omega .= scale * randn(sys.dims)
    end
    residuals!(sys)
    return nothing
end

function fieldinit!(sys::FieldSystem, fields::Dict; field_type::Symbol = :omega, interpolate::Bool = true)
    if field_type == :omega
        dict = sys.fields
    else
        dict = sys.density
    end

    for (mid, arr) in fields
        if haskey(sys.monomers, mid)
            if size(arr) == sys.dims
                dict[mid] = arr
            elseif interpolate
                itp = interpolate_grid(arr, sys.dims)
                dict[mid] = itp
            else
                error("Input field dimensions $(size(arr)) does not match system grid $(sys.dims).")
            end
        else
            @warn "Monomer with ID = $(mid) not present in system -- skipping."
        end
    end 

    # Update based on fields
    if field_type == :rho
        uniform_fields!(sys)
    end
    residuals!(sys)

    return nothing
end

function fieldinit!(sys::FieldSystem, fpath::AbstractString; load_cell::Bool = true, interpolate::Bool = true)
    data = loadfields(fpath)
    if load_cell
        cell = UnitCell(data[:dim], data[:crystal], data[:cell_params])
        add_cell!(sys, cell)
    end
    fieldinit!(sys, data[:fields]; field_type = data[:field_type], interpolate = interpolate)
    return nothing
end


#==============================================================================#
# Field IO
#==============================================================================#

"""
    savefields(sys, fpath; type = :field/:density)

Save fields from the system to a JSON file.
Argument `type` specifies whether density or omega fields are saved.
"""
function savefields(sys::FieldSystem, fpath::AbstractString; field_type::Symbol = :omega)  
    d = Dict()
    d[:dim] = ndims(sys.dims)
    d[:grid] = sys.dims
    d[:mids] = sort(collect(keys(sys.monomers)))

    d[:date] = Dates.format(Dates.now(), "Y-m-d H:M:S")
    d[:crystal] = sys.cell.crystal
    d[:cell_params] = sys.cell.params
    d[:field_type] = field_type

    fields = field_type == :rho ? sys.density : sys.fields
    d[:fields] = typeof(fields)()
    for (mid, arr) in fields
        d[:fields][mid] = arr
    end

    json_string = JSON.json(d, 4)#JSON.json(d, 4)
    open(fpath, "w") do f 
        write(f, json_string) 
    end

    return nothing
end

"""
    loadfields(fpath)

Load fields and system data from a specified JSON file. 
"""
function loadfields(fpath::AbstractString)
    raw = JSON.parsefile(fpath)

    d = Dict()
    d[:dim] = Int(raw["dim"])
    d[:grid] = Tuple(convert(Vector{Int}, raw["grid"]))

    d[:mids] = Tuple(convert(Vector{Int}, raw["mids"]))
    d[:crystal] = Symbol(raw["crystal"])
    d[:cell_params] = convert(Vector{Float64}, raw["cell_params"])
    d[:field_type] = Symbol(raw["field_type"])

    mids = d[:mids]
    grid = d[:grid]
    d[:fields] = Dict(mid => zeros(grid) for mid in mids)
    for mid in mids
        arr = d[:fields][mid]
        vecs = raw["fields"][string(mid)]
        for k = 1:grid[end]
            for j = 1:grid[2]
                arr[:,j,k] = vecs[k][j]
            end
        end
    end

    return d
end