module StraightStellarator

using DifferentialEquations
using Plots
using QuadGK
using StaticArrays

export μ0, Coil, CoilConfiguration, compute_B, compute_ψ
export compute_poincare, poincare_plot
export cross, cross!

include("Coil.jl")
include("utils.jl")
include("results.jl")

end
