"""
	struct Monomer

A chemical monomer type, identified by a unique monomer ID
that is used to access a field and density grid from within a system.
"""
struct Monomer
	id     :: Int
	size   :: Float64
	charge :: Float64
    name   :: String

	function Monomer(; id::Integer, size::Real = 1.0, charge::Real = 0.0, name::String = "")
		return new(id, size, charge, name)
	end
end

#==============================================================================#

function Base.Base.show(io::IO, mon::Monomer)
    if !isempty(mon.name)
        @printf(io, "Monomer(`%s`, mid = %d, size = %.2f)", mon.name, mon.id, mon.size)
    else
        @printf(io, "Monomer(mid = %d, size = %.2f)", mon.id, mon.size)
    end
end

Base.hash(mon::Monomer, h::UInt) = hash(mon.id, h)

Base.isequal(amon::Monomer, bmon::Monomer) = (hash(amon) == hash(bmon))