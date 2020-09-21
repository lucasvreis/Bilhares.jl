__precompile__()
module BilharesMarkov

export prob_estacionária, mat_estocástica, map_markov, média_temporal

using DynamicalBilliards
import DynamicalBilliards: normalvec, specular!, collision, extrapolate, cossin,
                           nocollision, ispinned, accuracy, check_condition, timeseries!
using StaticArrays, LinearAlgebra, Arpack, SparseArrays
using ProgressMeter
using PyPlot
const SV = SVector{2}
const SV3 = SVector{3}
const SM = SMatrix{2,2}

include("utils/utils.jl")
include("linearlypropagating.jl")
include("hyperbolic/hyperbolic.jl")
include("elliptic.jl")
include("analysis/markov.jl")
include("analysis/monte-carlo-integration.jl")

end
