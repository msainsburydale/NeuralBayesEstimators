# The first convolution layer transforms the input array, which is of
# size (29, 23), into an array of size (16, 16). Since this transformed array
# is the same size as the assumed square grid in the simulation studies, we may
# then set the remaining layers to neural network used in the simulation studies.
# Note that the formula to determine the appropriate kernel size is, assuming a
# kernal stride (how many pixels the kernel moves at each step) of 1, we have
# 			hₒ = hᵢ - hₖ + 1,
# 			wₒ = wᵢ - wₖ + 1,
# where hₒ is the height of output array, hᵢ is the height of the input array,
# hₖ is the height of the convolutional kernal, and similarly for wₒ, wᵢ, and wₖ,
# where w stands for width.

using Flux
using SpatialDeepSets: loadwithoutdict

path = joinpath(pwd(), "data/RedSea/regular/")
wᵢ = loadwithoutdict(path * "width.rda", "width")
hᵢ = loadwithoutdict(path * "height.rda", "height")
wₒ = 16
hₒ = 16
wₖ = wᵢ - wₒ + 1
hₖ = hᵢ - hₒ + 1

function architecture(p)

	qₜ = 128

	ψ = Chain(
		Conv((hₖ, wₖ), 1 => 16, relu),
		Conv((10, 10), 16 => 32,  relu),
		Conv((5, 5),  32 => 64,  relu),
		Conv((3, 3),  64 => qₜ, relu),
		Flux.flatten
		)

	ϕ = Chain(
		Dense(qₜ, 500, relu),
        Dense(500, p),
		finallayer
	)

	return ψ, ϕ
end
architecture() = architecture(p)

function finallayer(x)
	θ̂₁ = exp.(x[1:5, :])
	θ̂₂ = x[6, :]' # NB This is the idx for μ, which I've just hard-coded for simplicity.
	θ̂₃ = exp.(x[7:8, :])
	return vcat(θ̂₁, θ̂₂, θ̂₃)
end
