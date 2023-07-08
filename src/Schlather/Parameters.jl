using NeuralEstimators
using SpatialDeepSets
using Distances
using Distributions: Uniform

Ω = (
	ρ = Uniform(2.0, 10.0),
	ν = Uniform(0.5, 3.0)
)
parameter_names = String.(collect(keys(Ω)))
S = expandgrid(1:16, 1:16)
S = Float64.(S)
D = pairwise(Euclidean(), S, S, dims = 1)
ξ = (
	Ω = Ω, S = S, D = D, p = length(Ω),
	parameter_names = parameter_names,
	ρ_idx = findfirst(parameter_names .== "ρ"),
	ν_idx = findfirst(parameter_names .== "ν"),
	invtransform = exp
)
