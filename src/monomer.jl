"""
	struct Monomer

A chemical monomer type, identified by a unique monomer ID
that is used to access a field and density grid from within a system.
"""
struct Monomer
	id     :: Int
	size   :: Float64
	charge :: Float64

	function Monomer(; id::Integer, size::Real = 1.0, charge::Real = 0.0)
		return new(id, size, charge)
	end
end

#==============================================================================#

Base.show(io::IO, mon::Monomer) = @printf(io, "Monomer(mid = %d, size = %.2f)", mon.id, mon.size)

Base.hash(mon::Monomer, h::UInt) = hash(mon.charge, hash(mon.size, hash(mon.id, h)))

Base.:(==)(a::Monomer, b::Monomer) = (hash(a) == hash(b))

Base.isless(a::Monomer, b::Monomer) = (a.id < b.id)