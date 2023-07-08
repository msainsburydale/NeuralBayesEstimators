using NeuralEstimators
using SpatialDeepSets
using Distributions: Normal, mean

# prior distribution of θ
if prior == "diffuse"
	@info "Using a diffuse prior"
	Ω = Normal(0, 5)
elseif prior == "informative"
	@info "Using an informative prior"
	Ω = Normal(0, 0.5)
else
	error("Unrecognised prior")
end

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
