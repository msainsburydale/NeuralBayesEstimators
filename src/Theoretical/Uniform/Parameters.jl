using NeuralEstimators
using SpatialDeepSets
using Distributions: Uniform, Pareto, mean, quantile

# prior distribution for θ
Ω = Pareto(4, 1)

ξ = (
	Ω = Ω,
	parameter_names = ["θ"]
)

struct Parameters <: ParameterConfigurations
	θ
end

Parameters(K::Integer, ξ) = Parameters(rand(ξ.Ω, 1, K))

θ_scenarios = mean(ξ.Ω)
θ_scenarios = reshape([Float32(θ_scenarios)], 1, 1)
θ_scenarios = Parameters(θ_scenarios)
