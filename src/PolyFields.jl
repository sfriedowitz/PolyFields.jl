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
using Parameters

using Random
using StaticArrays
using Interpolations
using FFTW

#==============================================================================#
# Exports
#==============================================================================#

# Types
export AbstractSystem, AbstractSpecies, AbstractConstraint, AbstractPropagator, AbstractInteraction, AbstractFieldUpdater
export FieldGrid, FieldDims, PseudoSpectralSolver, FFTHolder, UnitCell
export Monomer, Point, Homopolymer, Diblock, Multiblock
export FloryInteraction, EdwardsInteraction
export Ensemble, Canonical, Grand, Compressibility
export FieldSystem, SCFTOptions, SCFTResults
export DomainUpdater, EulerUpdater, AndersonUpdater, MomentumUpdater

# Methods
export make_basis, basis_matrix, basis_dimensions
export add_monomer!, add_species!, add_cell!, add_interaction!, validate
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

# Cell and system
include("system/fft.jl")
include("system/mde.jl")
include("system/monomer.jl")
include("system/basis.jl")
include("system/cell.jl")
include("system/compress.jl")
include("system/system.jl")

# Species architectures
include("species/base.jl")
include("species/point.jl")
include("species/multiblock.jl")

# # Interaction classes
include("interactions/base.jl")
include("interactions/flory.jl")
include("interactions/edwards.jl")

# # Field updaters
include("updaters/base.jl")
include("updaters/euler.jl")
include("updaters/anderson.jl")
include("updaters/momentum.jl")
include("updaters/domain.jl")

include("fileio.jl")
include("simulate.jl")
# # include("viz.jl")

end