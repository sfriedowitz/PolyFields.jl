# PolyFields.jl

PolyFields is a minimal Julia implementation of polymer self-consistent field theory,
as described in the extensive monograph
"Fredrickson, G. (2005). The Equilibrium Theory of Inhomogeneous Polymers."
Core implementation details are heavily inspired by the Fortran package [PSCF](https://github.com/dmorse/pscf).

Currently, the package only has the capability of solving the self-consistent field theory (SCFT) equations,
corresponding to a mean-field solution of the underlying field-theoretic partition function.
Implementation of full field-theoretic simultions (FTS) via complex Langevin sampling is a later goal.

## Usage

The PolyFields package is amenable for use in an interactive workspace environment.
The basic structs required for a field-theoretic simulation can be created as follows:

```julia
using PolyFields

# Number of grid-points in x/y/z dimensions
dims = (64, 1, 1);

# Unit cell for 1D system
cell = UnitCell(1, :lamellar, 10.0); 

# Create monomers, which represent a distinct chemical type
amon = Monomer(; id = 1, vol = 1.0, charge = 0.0);
bmon = Monomer(; id = 2, vol = 1.0, charge = 0.0);

# Create a diblock copolymer species
# Chain structure: N = 100, fA = 0.5, bA = bB = 1.0, Ns = 100
chain = Diblock(amon, bmon, 100, 0.5, 1.0, 1.0, 100);

# Create a Flory-Huggins interaction
itx = FloryInteraction();
set_interaction!(itx, amon, bmon, 0.2)

# Create a system to hold all the pieces
sys = FieldSystem(dims, cell; monomers = [amon, bmon], ensemble = Canonical);
add_species!(sys, chain)
add_interaction!(sys, itx)

# Initialize the fields randomly with a scaling factor of 0.1
fieldinit!(sys; scale = 0.1)
```

## TODO

* Add plotting tools for system and field output
* GCE performance and stability
* Add FDDD2 crystal systems, more generic format for basis construction