using Flux

function architecture(p, qₛ = 0)

	qₜ = 256

	ψ = Chain(
		Conv((10, 10), 1 => 64,  relu),
		Conv((5, 5),  64 => 128,  relu),
		Conv((3, 3),  128 => qₜ, relu),
		Flux.flatten
		)

	ϕ = Chain(
		Dense(qₜ + qₛ, 500, relu),
		Dense(500, p)
	)

	return ψ, ϕ
end
