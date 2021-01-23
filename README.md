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
```

## URGENT TODO

* Use symbols instead of strings for parameter flags
* Change method names to more Julia-like naming conventions (single words)
* Constructors for Cell should be changed
* Fix constraints and application during field optimizer
* Clean up field optimizers and interaction structure, referencing to system, method signatures
* Revamp field file format -- structured data format
* Add plotting functionality for 1/2/3D

## Roadmap

* Variable cell implementation working and stable
* Revamp structure of interactions and field updaters -- add abstract methods and pass system instead of binding
* Implementation of LDA approach for Gaussian electrostatics
* Canonical vs GCE calculation working accurately
* Monomer reference volume scaling tested and accurate (invariant to v and N shifts)
* Expanded IO methods for seeding, file output, logging, repository of stable seeds
* Expanded field updater methods -- working conjugate gradient and line search methods
* Compile PSCF basis functions, incorporate symmetry
* Optional density smearing methods
* Expansion of interactions to true electrostatics, dielectric inhomogeneity (iPSCF analogue)
* Complex langevin field updating schemes
