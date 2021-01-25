"""
# PolyFields.jl
Welcome to PolyFields.jl!
"""
module PolyFields

#==============================================================================#
# Dependencies
#==============================================================================#

import Statistics: mean, std
import LinearAlgebra: norm, dot, det, mul!, inv!, lu!
import Dates

using Printf
using JSON

using Random
using StaticArrays
using Interpolations
using FFTW
using LineSearches

#==============================================================================#
# Exports
#==============================================================================#

# Types
export AbstractSystem, AbstractSpecies, AbstractConstraint, AbstractPropagator, AbstractInteraction, AbstractFieldUpdater
export FieldGrid, FieldDims, Cell, Cell1D, Cell2D, Cell3D, FFTHolder, PseudoSpectralSolver
export Monomer, Point, Homopolymer, Diblock, Multiblock
export FloryInteraction, EdwardsInteraction
export Ensemble, Canonical, Grand, Compressibility, Electroneutrality
export FieldSystem, SCFTOptions, SCFTResults
export DomainUpdater, EulerUpdater, AndersonUpdater, MomentumUpdater

# Methods
export add_monomer!, add_species!, add_interaction!, add_constraint!, validate
export density!, muphi!, potentials!, potential!, scfstress, scfstress!, uniform_fields!, residuals!
export free_energy, free_energy_bulk, energy, energy_bulk, set_interaction!
export step!, field_error, stress_error
export fieldinit!, savefields, loadfields, scft!

#==============================================================================#
# Load files
#==============================================================================#

# General
include("types.jl")
include("utils.jl")
include("fft.jl")
include("mde.jl")
include("cell.jl")
include("monomer.jl")

# System fields
include("system.jl")

# Species architectures
include("species/base.jl")
include("species/point.jl")
include("species/multiblock.jl")

# # Interaction classes
include("interactions/base.jl")
include("interactions/flory.jl")
include("interactions/edwards.jl")

# # Field updaters
include("constraints.jl")
include("updaters/base.jl")
# include("updaters/constraint.jl")
# include("updaters/domain.jl")
# include("updaters/euler.jl")
# include("updaters/anderson.jl")
# include("updaters/momentum.jl")
# # include("updaters/linesearch.jl")
# # include("updaters/conjgrad.jl")

include("fileio.jl")
#include("simulate.jl")
# include("plots.jl")

end