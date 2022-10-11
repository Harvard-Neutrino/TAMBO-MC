module Tambo

export Injector

using StaticArrays
using Rotations
using Random
using StatsBase: sample
using PyCall
using LinearAlgebra:norm
using Roots: find_zeros, find_zero
using JLD2
using Dierckx
using Distributions: Uniform
using ProgressBars

include("units.jl")
include("powerlaws.jl")
include("geometries.jl")
include("tracks.jl")
include("crosssections.jl")
include("particles.jl")
include("proposalinterface.jl")
include("inject.jl")

mutable struct Simulator
    
end

end # module
