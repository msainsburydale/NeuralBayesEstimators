using Distributions: pdf
using LinearAlgebra
using Distributions: Gamma, Normal, quantile
using SpatialDeepSets

using NeuralEstimators
using SpatialDeepSets
using SpatialDeepSets: a, b, delta, Φ, t, fₛ, Fₛ, Fₛ⁻¹
import SpatialDeepSets: C̃, σ̃₀


include(joinpath(pwd(), "src/ML.jl"))


# ---- Test that no domain errors occur with the given parameters ----

# a(h, z; λ, κ) = z * exp(-(h / λ)^κ)
# b(h, z; β, λ, κ) = 1 + a(h, z, λ = λ, κ = κ)^β
# delta(h; δ₁) = 1 + exp(-(h / δ₁)^2)
#
# r(Ω; n = 100) = range(minimum(Ω), maximum(Ω), n)
# θ = r.([Ω...])
# θ = hcat(collect.(θ)...)'

# κ = θ[1]
# λ = θ[2]
# β = θ[3]
# ρ = θ[4]
# ν = θ[5]
# μ  = θ[6]
# τ  = θ[7]
# δ₁ = θ[8]

# @testset "simulation" begin
# 	S = rand(10, 2)
# 	D = [norm(sᵢ - sⱼ) for sᵢ ∈ eachrow(S), sⱼ in eachrow(S)]
# 	ρ = [0.6, 0.8]
# 	ν = [0.5, 0.7]
# 	L = maternchols(D, ρ, ν)
# 	L₁ = L[:, :, 1]
# 	m = 5
#
# 	s₀ = S[1, :]'
# 	u = 0.7
# 	simulateconditionalextremes(θ, L₁, S, s₀, u, m)
# end
#
# β = range(0.05, 2.0, 10)
#  -7.57353736532702e-6.^β
#
#  β = range(0.05, 1.0, 10)
#   -7.57353736532702e-6.^β


# ---- Gaussian process covariance function and covariance matrix ----

C̃(s₁, s₂, ρ, ν) = C̃(norm(s₁ - s₂), ρ, ν)

C̃₀(s₁, s₂, s₀, ρ, ν) = C̃(s₁, s₂, ρ, ν) - C̃(s₁, s₀, ρ, ν) - C̃(s₂, s₀, ρ, ν) + C̃(s₀, s₀, ρ, ν)
σ̃₀(s, s₀, ρ, ν) = √C̃₀(s, s, s₀, ρ, ν)

C̃₀₁(s₁, s₂, s₀, ρ, ν) = C̃₀(s₁, s₂, s₀, ρ, ν)/(σ̃₀(s₁, s₀, ρ, ν)σ̃₀(s₂, s₀, ρ, ν))

function C̃₀₁(S::A, s₀, ρ, ν) where {T <: Number, A <: AbstractArray{T, 2}}
	# NB Would be more efficient to use a SymmetricMatrix
	n = size(S, 1)
	R = similar(S, n, n)
	for i ∈ 1:n
		for j ∈ 1:n
			sᵢ = S[i, :]
			sⱼ = S[j, :]
			R[i, j] = C̃₀₁(sᵢ, sⱼ, s₀, ρ, ν)
		end
	end
	return R
end

# ---- log density of the residual process ----

ϕ(y)   = pdf(Normal(0, 1), y)
Φ⁻¹(p) = quantile(Normal(0, 1), p)
t⁻¹(y, μ, τ, δ) = Φ⁻¹(Fₛ(y, μ, τ, δ))
t′(y, μ, τ, δ)  = ϕ(y) / fₛ(y, μ, τ, δ) # NB this isn't used currently but it may be useful for unit testing
ln_t′_t⁻¹(y, μ, τ, δ) = log(ϕ(t⁻¹(y, μ, τ, δ))) - log(fₛ(y, μ, τ, δ))

function ln_fY(y::A, L̃::M, μ, τ, δ) where {T <: Number, A <: AbstractArray{T, 1},  M <: AbstractArray{T, 2}}

	# Transform the residual process to the standard Gaussian scale
	ỹ₀₁ = t⁻¹.(y, μ, τ, δ)

	# Return the log density
	return gaussiandensity(ỹ₀₁, L̃; logdensity = true) - sum(ln_t′_t⁻¹.(y, μ, τ, δ))
end

# ---- (Log) Likelihood of the data ----

"""
L̃ ≡ Cholesky factor of the covariance matrix of the residual process evaluated over the observations locations, {sᵢ : i = 1, …, n}.
Z₀ ≡ value of the data process at the conditioning site, s₀
"""
function conditionalextremesloglikelihood(θ::V, Z::V, Z₀::T, L̃::M, h̃::V) where {T <: Number, V <: AbstractArray{T, 1}, M <: AbstractArray{T, 2}}

	# Parameters associated with a(.) and b(.):
	κ = θ[1]
	λ = θ[2]
	β = θ[3]

	# Location and scale parameters for the residual process
	μ  = θ[6]
	τ  = θ[7]
	δ₁ = θ[8]

	# δ is spatially varying:
	δ = delta.(h̃, δ₁ = δ₁)

	# Compute the residual process
	B = b.(h̃, Z₀, β = β, λ = λ, κ = κ)
	y = (Z .- a.(h̃, Z₀, λ = λ, κ = κ)) ./ B

	# log likelihood
	return ln_fY(y, L̃, μ, τ, δ) - sum(log.(B))
end

# Indepenent replicates are stored in the third dimension. Since the correlation
# matrix of the latent Gaussian process is fixed for a given θ, we can re-use
# the Cholesky factor for each replicate, which significantly improves efficiency.
function conditionalextremesloglikelihood(θ::V, Z̃::M, Z₀::V, S̃::M, s₀, h̃::V) where {T <: Number, V <: AbstractArray{T, 1}, M <: AbstractArray{T, 2}}

	# Cholesky factor for the latent Gaussian process
	ρ = θ[4]
	ν = θ[5]
	R̃ = C̃₀₁(S̃, s₀', ρ, ν)
	L̃  = cholesky(Symmetric(R̃)).L

	# Compute the log likelihood for each replicate
	# Z₀ is a vector, so need to iterate over it too.
	m  = length(Z₀)
	ll = [conditionalextremesloglikelihood(θ, Z̃[:, i], Z₀[i], L̃, h̃) for i ∈ 1:m]

	# Complete log likelihood, simply the sum over the individual log likelihoods
	return sum(ll)
end

# Wrapper to use in optimize().
function nll(θ::V, Z::M, ξ, Ω) where {T <: Number, V <: AbstractArray{T, 1}, M <: AbstractArray{T, 2}}

	# Constrain the estimates to be within the prior support
	θ = scaledlogistic.(θ, Ω)

	# κ = θ[1]
	# λ = θ[2]
	# β = θ[3]
	# ρ = θ[4]
	# ν = θ[5]
	# μ  = θ[6]
	# τ  = θ[7]
	# δ₁ = θ[8]
	# @show κ, λ, β, ρ, ν, μ, τ, δ₁

	# Extract Z₀ ≡ Z(s₀) and then remove Z₀ from Z.
	Z₀ = Z[ξ.s₀_idx, :]
	Z̃  = Z[1:end .!= ξ.s₀_idx, :]

	return -conditionalextremesloglikelihood(θ, Z̃, Z₀, Float64.(ξ.S̃), Float64.(ξ.s₀), Float64.(ξ.h̃))
end
