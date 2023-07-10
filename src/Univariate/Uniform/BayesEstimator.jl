function MLE(Z) where {T <: Number, N <: Int, A <: AbstractArray{T, N}, V <: AbstractVector{A}}
    maximum.(Z)'
end

MLE(Z, ξ) = MLE(Z)

# ---- Bayes estimator under the conjugate prior ----

using Distributions: median, shape, scale
function BayesEstimator(Z, ξ) where {T <: Number, N <: Int, A <: AbstractArray{T, N}, V <: AbstractVector{A}}
    median.(posterior.(Z, Ref(ξ)))'
end

function posterior(Z, ξ)
    α = shape(ξ.Ω)
    β = scale(ξ.Ω)
    Pareto(α̃(α, Z), β̃(β, Z))
end

α̃(α, Z) = α + length(Z)
β̃(β, Z) = maximum([Z..., β])


# ---- One-at-a-time estimator under the conjugate prior ----

function OAATEstimator(Z, ξ) where {T <: Number, N <: Int, A <: AbstractArray{T, N}, V <: AbstractVector{A}}
	θ̂ = map(Z) do z
		mean(median.(posterior.(z, Ref(ξ))))
	end
	θ̂'
end

# ---- MAP estimator under the conjugate prior ----

using Distributions: mode
function MAPEstimator(Z, ξ) where {T <: Number, N <: Int, A <: AbstractArray{T, N}, V <: AbstractVector{A}}
    mode.(posterior.(Z, Ref(ξ)))'
end
