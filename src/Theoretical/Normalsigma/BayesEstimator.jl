using Distributions: InverseGamma, median, shape, scale
residualsumsquares(Z, μ) = sum((Z .- μ).^2)

# ---- Bayes estimator under the conjugate prior ----

function BayesEstimator(Z, ξ) where {T <: Number, N <: Int, A <: AbstractArray{T, N}, V <: AbstractVector{A}}
    median.(posterior.(Z, Ref(ξ)))'
end

function posterior(Z, ξ)
    μ = ξ.μ
    Ω = ξ.Ω
    α = shape(Ω)
    β = scale(Ω)
    m = length(Z)

    InverseGamma(α̃(α, m), β̃(β, Z, μ))
end

α̃(α, m) = α + m/2
β̃(β, Z, μ) = β + residualsumsquares(Z, μ)/2

# ---- MLE ----

"""
See [here](https://math.stackexchange.com/a/2187053) for the derivation of MLE.
"""
function MLE(Z, ξ) where {T <: Number, N <: Int, A <: AbstractArray{T, N}, V <: AbstractVector{A}}

    μ = ξ.μ

    θ̂ = [residualsumsquares(z, μ) / length(z) for z ∈ Z]

    reshape(θ̂, 1, :) # return as p x K matrix
end
