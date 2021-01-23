"""
# PolyFields.jl
Welcome to PolyFields.jl!

PolyFields is a package...

## REPL help

## Documentation

"""
module PolyFields

#==============================================================================#
# Dependencies
#==============================================================================#

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
# Exports
#==============================================================================#

# Types
export AbstractSpecies, AbstractConstraint, AbstractInteraction, AbstractFieldUpdater
export FieldGrid, Cell, Cell1D, Cell2D, Cell3D, FFTBuddy, PseudoSpectralSolver
export Monomer, Point, Homopolymer, Diblock, Multiblock
export ChiInteraction, GaussianEdwards
export Ensemble, Canonical, Grand, Compressibility, Electroneutrality
export FieldSystem, SCFTOptions, SCFTResults
export DomainUpdater, EulerUpdater, AndersonUpdater, MomentumUpdater

# Methods
export add_mde_solver!, add_monomer!, add_species!, add_interaction!, validate
export density!, muphi!, potentials!, potential!, update_eta!, stress!, set_omega_uniform, residuals!
export free_energy, free_energy_bulk, energy, energy_bulk, set_interaction!
export step!, field_error, stress_error
export scft!, fieldinit!, fieldsave, fieldload!

#==============================================================================#
# Load files
#==============================================================================#

# General
include("types.jl")
include("utils.jl")

include("cell.jl")
# include("fourier.jl")
# include("mde.jl")
# include("system.jl")

# # Species architectures
# include("species/base.jl")
# include("species/monomer.jl")
# include("species/point.jl")
# include("species/multiblock.jl")

# # Interaction classes
# include("interactions/base.jl")
# include("interactions/chi.jl")
# include("interactions/gaussian.jl")

# # Field updaters
# include("updaters/base.jl")
# include("updaters/constraint.jl")
# include("updaters/domain.jl")
# include("updaters/euler.jl")
# include("updaters/anderson.jl")
# include("updaters/momentum.jl")
# # include("updaters/linesearch.jl")
# # include("updaters/conjgrad.jl")

# include("simulate.jl")
# include("file_io.jl")
# include("plots.jl")

end