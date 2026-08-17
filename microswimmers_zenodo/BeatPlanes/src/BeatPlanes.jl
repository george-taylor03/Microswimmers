module BeatPlanes

using MicroSwimmers
using JLD2
using Parameters
using Statistics
using OrdinaryDiffEq
using DiffEqCallbacks
using StaticArrays
using LinearAlgebra

include("simulations.jl")
include("chlamy.jl")
include("filter_feeder.jl")
include("excavate_2_chlamy.jl")
include("exports.jl")

end # module BeatPlanes
