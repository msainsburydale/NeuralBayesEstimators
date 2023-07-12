using Flux
function architecture(n::Integer, p::Integer)

	qₜ = 256

	ψ = Chain(
		Dense(n, 300, relu),
		Dense(300, qₜ, relu)
	)

	ϕ = Chain(
		Dense(qₜ, 50, relu),
		Dense(50, p),
		x -> dropdims(x, dims = 2), # Drop the "channel" dimension, which is a singelton dimension; NB: Would be better to change the data simulation to not include a channel dimension
		finallayer
	)

	return ψ, ϕ
end
architecture() = architecture(n, p)

function finallayer(x)
	θ̂₁ = exp.(x[1:5, :])
	θ̂₂ = x[6, :]' # NB This is the idx for μ, which I've just hard-coded for simplicity.
	θ̂₃ = exp.(x[7:8, :])
	return vcat(θ̂₁, θ̂₂, θ̂₃)
end
