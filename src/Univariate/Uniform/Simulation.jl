import NeuralEstimators: simulate
using Distributions: Uniform

function simulate(params::Parameters, m::R) where {R <: AbstractRange{I}} where I <: Integer
	K = size(params, 2)
	m̃ = rand(m, K)
	θ = params.θ

	return [rand(Uniform(0, θ[k]), 1, 1, m̃[k]) for k ∈ eachindex(m̃)]
end
simulate(parameters::Parameters, m::Integer) = simulate(parameters, range(m, m))
