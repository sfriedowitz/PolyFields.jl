"""
# PolyFields.jl
Welcome to PolyFields.jl!

PolyFields is a package...

## REPL help

## Documentation

"""
module PolyFields

#==============================================================================#
# Exports
#==============================================================================#

# Types
export AbstractSpecies, AbstractConstraint, AbstractInteraction, AbstractFieldUpdater
export NPWGrid, Cell, Cell1D, Cell2D, Cell3D, FFTBuddy, PseudoSpectralSolver
export Monomer, Point, Homopolymer, Diblock, Multiblock
export FloryInteraction, GaussianEdwards
export Ensemble, Canonical, Grand, Compressibility, Electroneutrality, FieldSystem, SCFTOptions
export DomainUpdater, EulerUpdater, AndersonUpdater, MomentumUpdater

# Methods
export add_mde_solver!, add_monomer!, add_species!, add_interaction!, validate
export update_density!, update_mu_phi!, update_potentials!, update_eta!, update_stress!, set_omega_uniform!, make_residuals!
export free_energy, free_energy_bulk, energy, energy_bulk, add_potential!, set_interaction!
export update_state!, field_error, stress_error
export scft!, init_fields!, save_fields, load_fields

#==============================================================================#
# Dependencies
#==============================================================================#

import Base: show, hash, isequal
import Statistics: mean, std
import LinearAlgebra: norm, dot, det, mul!, inv!, lu!
import Dates

using Parameters
using Printf
using JLD
using Random
using StaticArrays
using Interpolations
using FFTW
using LineSearches

#==============================================================================#
# Load files
#==============================================================================#

# General
include("types.jl")
include("utils.jl")

include("cell.jl")
include("monomer.jl")
include("fourier.jl")
include("mde.jl")
include("system.jl")
include("constraint.jl")

# Species architectures
include("species/abstract.jl")
include("species/point.jl")
include("species/multiblock.jl")

# Interaction classes
include("interactions/abstract.jl")
include("interactions/flory.jl")
include("interactions/gaussian.jl")

# Field updaters
include("updaters/abstract.jl")
include("updaters/domain.jl")
include("updaters/euler.jl")
include("updaters/anderson.jl")
include("updaters/momentum.jl")
# include("updaters/linesearch.jl")
# include("updaters/conjgrad.jl")

include("simulate.jl")
include("file_io.jl")
include("plots.jl")

end