module StraightStellarator

using QuadGK
using DifferentialEquations
using Plots

export μ0, Coil, CoilConfiguration, compute_B, compute_poincare

include("Coil.jl")
include("utils.jl")
include("results.jl")

end
