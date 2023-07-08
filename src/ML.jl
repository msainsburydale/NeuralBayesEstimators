using LinearAlgebra
using NeuralEstimators
using SpatialDeepSets
using Optim
using Folds # parallel version of map()
using Flux: flatten
using RecursiveArrayTools

# ---- Likelihood estimators ----

function likelihoodestimator(Z::M, θ₀::V, ξ, Ω) where {T <: Number, V <: AbstractArray{T, 1}, M <: AbstractArray{T, 2}}

	# Closure that will be minimised:
	loss(θ) = nll(θ, Z, ξ, Ω)

	# Estimate the parameters
	θ̂ = optimize(loss, θ₀, NelderMead()) |> Optim.minimizer

	# During optimisation, we constrained the parameters using the scaled-logistic
	# function; here, we convert to the orginal scale
	θ̂ = scaledlogistic.(θ̂, Ω)

	return θ̂
end


function likelihoodestimator(Z::V, ξ) where {T <: Number, M <: AbstractArray{T, 2}, V <: AbstractVector{M}}

	# intitialise the ML estimates to the true parameters. Since we logistic-transform
	# the parameters during optimisation to force the estimates to be within the
	# prior support, here we provide the logit-transformed values.
	Ω = ξ.Ω
	Ω = [Ω...] # convert to array since broadcasting over dictionaries and NamedTuples is reserved
	θ₀ = scaledlogit.(ξ.θ₀, Ω)

	# Inverse of variance-stabilising transformation:
	Z = broadcast.(ξ.invtransform, Z)

	# Convert to Float64 so that Cholesky factorisation doesn't throw positive
	# definite error due to rounding.
	# (When you want to broadcast broadcast, then broadcast broadcast)
	Z  = broadcast.(Float64, Z)
	θ₀ = Float64.(θ₀)

	# Z may be repeated: Repeat θ₀ accordingly.
	θ₀ = repeatθ₀(θ₀, Z)

	# Number of parameter configurations to estimate
	K = size(θ₀, 2)

	# Convert from matrix to vector of vectors
	θ₀ = [θ₀[:, k] for k ∈ 1:K]

	# Optimise
	θ̂ = Folds.map(Z, θ₀) do Zₖ, θ₀ₖ
		 likelihoodestimator(Zₖ, θ₀ₖ, ξ, Ω)
	end

	# Convert to matrix
	θ̂ = convert(Base.typename(M).wrapper, VectorOfArray(θ̂))

	return θ̂
end


function likelihoodestimator(Z::V, ξ) where {T <: Number, A <: AbstractArray{T, 4}, V <: AbstractVector{A}}

	# Compress the data from a 4-dimensional array to a 2-dimensional array
	Z = flatten.(Z)

	return likelihoodestimator(Z, ξ)
end


# helper functions
function repeatθ₀(θ₀, Z)

	K = size(θ₀, 2)
	m  = length(Z)
	if m != K
		if (m ÷ K) == (m / K)
			θ₀ = repeat(θ₀, outer = (1, m ÷ K))
		else
			error("The length of the data vector, m, and the number of parameter configurations, K, do not match; further, m is not a multiple of K, so we cannot replicate θ to match m.")
		end
	end

	return θ₀
end
