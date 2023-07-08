using NeuralEstimators
using SpatialDeepSets
using Distributions: InverseGamma, mean

# prior distribution of θ
if prior == "diffuse"
	@info "Using a diffuse prior"
	Ω = InverseGamma(0.7, 0.7)
elseif prior == "informative"
	@info "Using an informative prior"
	Ω = InverseGamma(2, 2)
else
	error("Unrecognised prior")
end

ξ = (
	Ω = Ω,
	μ = 0,
	parameter_names = ["θ"]
)

struct Parameters <: ParameterConfigurations
	θ
end

Parameters(K::Integer, ξ) = Parameters(rand(ξ.Ω, 1, K))

θ_scenarios = mean(ξ.Ω)
θ_scenarios = reshape([Float32(θ_scenarios)], 1, 1)
θ_scenarios = Parameters(θ_scenarios)
