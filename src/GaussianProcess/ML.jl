using LinearAlgebra
using NeuralEstimators
using SpatialDeepSets
include(joinpath(pwd(), "src", "ML.jl"))


function covariancematrix(D; σₑ, ρ, ν)
    Σ = matern.(D, ρ, ν)
	Σ[diagind(Σ)] .+= σₑ^2
    return Σ
end

"""
Negative log-likelihood function to be minimised using Optim. If length(θ) > 2,
the smoothness parameter, ν, is estimated; otherwise, it is fixed to 1.
"""
function nll(θ, Z, ξ, Ω)

	# Constrain the estimates to be within the prior support
	θ = scaledlogistic.(θ, Ω)

	# Convert distance matrix to FLoat64 to avoid positive definite errors due to rounding
	D = Float64.(ξ.D)

	# Construct the covariance matrix with the current parameters
	ν = length(θ) > 2 ? θ[3] : 1
	Σ = covariancematrix(D, σₑ = θ[1], ρ = θ[2], ν = ν)

	return -gaussiandensity(Z, Σ; logdensity = true)
end
