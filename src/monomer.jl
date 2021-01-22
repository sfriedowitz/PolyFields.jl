"""
	struct Monomer

A struct representing a chemical monomer type.
A `Monomer` is identified by its unique monomer `id`, 
which is used to access fields and density grids within a system.

### Constructors
```julia
Monomer(; id, size = 1.0, charge = 0.0, name = "")
```
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

function show(io::IO, mon::Monomer)
    if !isempty(mon.name)
        @printf(io, "Monomer(`%s`, mid = %d, size = %.2f)", mon.name, mon.id, mon.size)
    else
        @printf(io, "Monomer(mid = %d, size = %.2f)", mon.id, mon.size)
    end
end

hash(mon::Monomer, h::UInt) = hash(mon.id, h)

isequal(amon::Monomer, bmon::Monomer) = (hash(amon) == hash(bmon))