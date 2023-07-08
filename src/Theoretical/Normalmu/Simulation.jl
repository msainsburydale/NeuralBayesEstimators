import NeuralEstimators: simulate
using Distributions: Normal

function simulate(params::Parameters, m::R) where {R <: AbstractRange{I}} where I <: Integer
	K = size(params, 2)
	m̃ = rand(m, K)
	θ = vec(params.θ)

	return [rand(Normal(θ[k], ξ.σ), 1, 1, m̃[k]) for k ∈ eachindex(m̃)]
end
simulate(parameters::Parameters, m::Integer) = simulate(parameters, range(m, m))
