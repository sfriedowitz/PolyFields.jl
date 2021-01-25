"""
	struct Monomer

A chemical monomer type, identified by a unique monomer ID
that accesses a field and density grid from the system.
"""
struct Monomer
	id     :: Int
	vol    :: Float64
	charge :: Float64

	function Monomer(; id::Integer, vol::Real = 1.0, charge::Real = 0.0)
		return new(id, vol, charge)
	end
end

#==============================================================================#

Base.show(io::IO, mon::Monomer) = @printf(io, "Monomer(mid = %d, vol = %.2f, charge = %.2f)", mon.id, mon.vol, mon.charge)

Base.hash(mon::Monomer, h::UInt) = hash(mon.charge, hash(mon.vol, hash(mon.id, h)))

Base.:(==)(a::Monomer, b::Monomer) = (hash(a) == hash(b))

Base.isless(a::Monomer, b::Monomer) = (a.id < b.id)