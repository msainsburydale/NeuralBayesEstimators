using Flux

function architecture(n, qₜ, qₛ, p; w = 128, d = 1, activation = relu)
	ψ = Chain(
		Dense(n, w, activation),
		[Dense(w, w, activation) for _ in 1:d]...,
		Dense(w, qₜ, activation)
	)
	ϕ = Chain(
		Dense(qₜ + qₛ, w, activation),
		[Dense(w, w, activation) for _ in 1:d]...,
		Dense(w, p),
		Flux.flatten
	)
	return ψ, ϕ
end
