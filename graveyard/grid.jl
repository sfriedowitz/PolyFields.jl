"""
	struct FieldGrid{T <: Number}

Minimal wrapper around a 3-dimensional array for storage of discretized r-space or k-space fields.
"""
struct FieldGrid{T <: Number}
	ndims  :: Int
	active :: Vector{Int}
	data   :: Array{T,3}
end

function FieldGrid(data::Array{T,3}) where {T}
	active = findall(x -> x > 1, size(data))
	ndims = length(active)
	return FieldGrid(ndims, active, data)
end

FieldGrid(::Type{T}, size::NTuple{3,<:Integer}) where {T} = FieldGrid(zeros(T, size))

FieldGrid(::Type{T}, n1::TI, n2::TI, n3::TI) where {T, TI <: Integer} = FieldGrid(T, (n1, n2, n3))

FieldGrid(::Type{T}, n1::TI, n2::TI) where {T, TI <: Integer} = FieldGrid(T, (n1, n2, 1))

FieldGrid(::Type{T}, n1::TI) where {T, TI <: Integer} = FieldGrid(T, (n1, 1, 1))

#==============================================================================#

Base.show(io::IO, grid::FieldGrid{T}) where {T} = @printf(io, "FieldGrid{%s}(size = %s)", T, size(grid))

Base.size(grid::FieldGrid) = size(grid.data)

Base.axes(grid::FieldGrid) = axes(grid.data)

Base.copy(grid::FieldGrid) = FieldGrid(copy(grid.data))

Base.copy!(dest::FieldGrid{T}, src::FieldGrid{T}) where {T} = copy!(dest.data, src.data)

Base.:(+)(g1::FieldGrid, g2::FieldGrid) = FieldGrid(g1.data .+ g2.data)
Base.:(*)(g1::FieldGrid, g2::FieldGrid) = FieldGrid(g1.data .* g2.data)
Base.:(*)(grid::FieldGrid, c::Number) = FieldGrid(c*grid.data)
Base.:(*)(c::Number, grid::FieldGrid) = FieldGrid(c*grid.data)

function addto!(dest::FieldGrid{T}, src::FieldGrid{T}) where {T}
	@assert size(src) == size(dest)
	dest.data .+= src.data
	return nothing
end