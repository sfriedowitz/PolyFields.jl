#==============================================================================#
# Abstract types
#==============================================================================#

"""
    abstract type AbstractSpecies

Supertype for all molecular species in a field-theoretic simulation.
"""
abstract type AbstractSpecies end

"""
    abstract type AbstractConstraint

Supertype for all physical constraints applied to a system.
"""
abstract type AbstractConstraint end

"""
    abstract type AbstractInteraction

Supertype for all excess interactions in a field-theoretic simulation.
"""
abstract type AbstractInteraction end

"""
    abstract type AbstractFieldUpdater

Supertype for iteration algorithms to update fields in a field-theoretic simulation.
"""
abstract type AbstractFieldUpdater end

#==============================================================================#
# Type definitions
#==============================================================================#

# 3D array grid alias
const PWGrid{T} = Array{T, 3}

# Nullable union holder
const Option{T} = Union{Nothing,T}

# Available methods from LineSearches.jl package
const LSUnion = Union{MoreThuente, BackTracking, HagerZhang, StrongWolfe}

# 3D vectors and matrices
const SVec3D{T} = SVector{3,T}
const MVec3D{T} = MVector{3,T}
const Vec3D{T} = Union{SVec3D{T}, MVec3D{T}}

const SMat3D{T} = SMatrix{3,3,T,9}
const MMat3D{T} = MMatrix{3,3,T,9}
const Mat3D{T} = Union{SMat3D{T}, MMat3D{T}}