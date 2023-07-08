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


# # finallayer should be function that acts on a p x K matrix
# function architecture(p, qₛ, finallayer)
#
# 	qₜ = 256
#
# 	ψ = Chain(
# 		Conv((10, 10), 1 => 64,  relu),
# 		Conv((5, 5),  64 => 128,  relu),
# 		Conv((3, 3),  128 => qₜ, relu),
# 		Flux.flatten
# 		)
#
# 	ϕ = Chain(
# 		Dense(qₜ + qₛ, 500, relu),
#         Dense(500, p),
# 		finallayer
# 	)
#
# 	return ψ, ϕ
# end


# Z = rand(256, 20)
# ψ, ϕ = architecture(8, 0, finallayer)
# ϕ(Z)

# if model == "ConditionalExtremes"
# 	function finallayer(x)
# 		θ̂₁ = exp.(x[1:5, :])
# 		θ̂₂ = x[6, :]'  # NB This is the idx for μ, which I've just hard-coded for simplicity.
# 		θ̂₃ = exp.(x[7:8, :])
# 		return vcat(θ̂₁, θ̂₂, θ̂₃)
# 	end
# else
# 	finallayer = x -> exp.(x)
# end
