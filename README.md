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
The basic objects required for a field-theoretic simulation can be created as follows.

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

Step (0     /  1500): ferr = 2.714097e+00, fhelm = 9.11369e-02
Updating cell dimensions: serr = 4.44357e-03
Step (10    /  1500): ferr = 9.052774e-01, fhelm = 8.30728e-02
Updating cell dimensions: serr = 7.68904e-04
Step (20    /  1500): ferr = 5.502302e-01, fhelm = 7.86577e-02
Updating cell dimensions: serr = 5.13565e-04
Step (30    /  1500): ferr = 9.215338e-04, fhelm = 7.76757e-02
Updating cell dimensions: serr = 3.56919e-04
Step (40    /  1500): ferr = 1.334250e-01, fhelm = 7.82135e-02
Updating cell dimensions: serr = 2.44857e-04
Step (50    /  1500): ferr = 1.434354e-01, fhelm = 7.94343e-02
Updating cell dimensions: serr = 1.44504e-04
Step (60    /  1500): ferr = 2.511615e-01, fhelm = 8.02967e-02
Updating cell dimensions: serr = 6.87233e-05
Step (70    /  1500): ferr = 3.708123e-01, fhelm = 8.12030e-02
Updating cell dimensions: serr = 2.33257e-05
Step (80    /  1500): ferr = 4.503784e-01, fhelm = 8.18496e-02
Updating cell dimensions: serr = 6.33112e-06
Step (90    /  1500): ferr = 4.630755e-01, fhelm = 8.20975e-02
Step (100   /  1500): ferr = 1.470985e-02, fhelm = 8.19481e-02
Updating cell dimensions: serr = 3.06109e-05
Step (110   /  1500): ferr = 1.179689e-01, fhelm = 8.14978e-02
Updating cell dimensions: serr = 1.00249e-05
Step (120   /  1500): ferr = 1.180627e-01, fhelm = 8.13219e-02
Step (130   /  1500): ferr = 1.019987e-05, fhelm = 8.12514e-02
Updating cell dimensions: serr = 4.60256e-06
Step (140   /  1500): ferr = 1.755261e-02, fhelm = 8.13044e-02
Converged @ step = 148: ferr = 6.588308e-06, fhelm = 8.13189e-02
Converged cell parameters: [16.65150848436396]

SCFTResults
  time: Float64 0.1011729
  nsteps: Int64 148
  converged: Bool true
  ferr: Float64 3.1223576732052154e-6
  serr: Float64 4.610473878519985e-6
  fhelm: Float64 0.08132056917007947
```

## TODO

* Add plotting tools for system and field output
* GCE performance and stability
* Add FDDD2 crystal systems, more generic format for basis construction