using NeuralEstimators
using SpatialDeepSets
using Distributions: Normal, mean

# prior distribution for θ
Ω = Normal(0, 0.5)

ξ = (
	Ω = Ω,
	σ = 1, # known standard deviation of the data
	parameter_names = ["θ"]
)

struct Parameters <: ParameterConfigurations
	θ
end

Parameters(K::Integer, ξ) = Parameters(rand(ξ.Ω, 1, K))

θ_scenarios = mean(ξ.Ω)
θ_scenarios = reshape([Float32(θ_scenarios)], 1, 1)
θ_scenarios = Parameters(θ_scenarios)
