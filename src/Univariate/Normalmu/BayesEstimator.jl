using Statistics: mean
using Distributions: Normal, mean, median, std

# ---- Bayes estimator under the conjugate prior ----

function BayesEstimator(Z, ξ) where {T <: Number, N <: Int, A <: AbstractArray{T, N}, V <: AbstractVector{A}}
    median.(posterior.(Z, Ref(ξ)))'
end

function posterior(Z, ξ)
    μ₀  = mean(ξ.Ω)
    σ₀  = std(ξ.Ω)
    σ₀² = σ₀^2
    σ²  = ξ.σ^2
    Normal(μ̃(μ₀, σ₀², σ², Z), σ̃(σ₀², σ², Z))
end

μ̃(μ₀, σ₀², σ², Z) = (1/σ₀² + length(Z)/σ²)^-1 * (μ₀/σ₀² + sum(Z)/σ²)
σ̃(σ₀², σ², Z) = sqrt((1/σ₀² + length(Z)/σ²)^-1)

# ---- MLE ----

"""
The MLE is simply the sample maximum, which does not require `ξ`.
"""
function MLE(Z) where {T <: Number, N <: Int, A <: AbstractArray{T, N}, V <: AbstractVector{A}}
    mean.(Z)'
end
MLE(Z, ξ) = MLE(Z)
