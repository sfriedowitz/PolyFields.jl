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
The basic objects required for a field-theoretic simulation can be created as follows:

```julia
using PolyFields

# Number of grid-points in x/y/z dimensions
dims = (64, 1, 1);

# Unit cell for 1D system
cell = UnitCell(1, :lamellar, 10.0); 

# Create monomers, which represent a distinct chemical type
amon = Monomer(; id = 1, vol = 1.0, charge = 0.0);
bmon = Monomer(; id = 2, vol = 1.0, charge = 0.0);

# Create a diblock copolymer species with: N = 100, fA = 0.5, bA = bB = 1.0, Ns = 100
chain = Diblock(amon, bmon, 100, 0.5, 1.0, 1.0, 100);

# Create a Flory-Huggins interaction
itx = FloryInteraction();
set_interaction!(itx, amon, bmon, 0.2);

# Create a system to hold all the pieces
sys = FieldSystem(dims, cell; monomers = [amon, bmon], ensemble = Canonical);
add_species!(sys, chain);
add_interaction!(sys, itx);

# Initialize the fields randomly with a scaling factor of 0.1
fieldinit!(sys; scale = 0.1)
```

Then, after creating a system and adding the appropriate components, we can run an SCFT calculation.
This is done by specifying an `AbstractFieldUpdater` to relax the field configurations,
and optionally the cell configurations as well.

```julia
# Runtime options for SCFT
opts = SCFTOptions(nsteps = 1500, nsout = 100, ftol = 1e-4, stol = 1e-4);

# Anderson mixing field updater, with a history size of 10 previous configurations
updater = AndersonUpdater(; nhist = 10);

# Variable cell implementation. Update cell every 5 SCFT steps, and perform 20 iterations
domain = DomainUpdater(; nskip = 5, nper = 20);

# Run the SCFT simulation
scft!(sys, updater; domain = domain, opts = opts)
Step (0     /  1500): ferr = 2.675499e+00, fhelm = 8.46080e-02
Updating cell dimensions: serr = 3.85896e-03
Step (10    /  1500): ferr = 3.178666e-01, fhelm = 8.11966e-02
Updating cell dimensions: serr = 2.63740e-04
Step (20    /  1500): ferr = 1.383014e-02, fhelm = 7.98100e-02
Updating cell dimensions: serr = 1.35652e-04
Step (30    /  1500): ferr = 3.056864e-02, fhelm = 8.07546e-02
Updating cell dimensions: serr = 3.80441e-05
Step (40    /  1500): ferr = 4.651606e-02, fhelm = 8.09580e-02
Updating cell dimensions: serr = 1.23038e-05
Step (50    /  1500): ferr = 5.035081e-02, fhelm = 8.12200e-02
Converged @ step = 55: ferr = 7.174001e-06, fhelm = 8.12299e-02
Converged cell parameters: [16.516666566336966]
```

## TODO

* Add plotting tools for system and field output
* GCE performance and stability
* Add FDDD2 crystal systems, more generic format for basis construction