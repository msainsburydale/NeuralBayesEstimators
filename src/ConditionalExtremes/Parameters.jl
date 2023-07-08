using Distances
using Distributions: Uniform
using NeuralEstimators
using SpatialDeepSets
using LinearAlgebra

standardLaplacequantile(p) = p < 0.5 ? log(2p) : -log(2 - 2p)

Ω = (
	# parameters associated with a(.) and b(.)
	κ = Uniform(1.0, 2.0),
	λ = Uniform(2.0, 5.0),
	β = Uniform(0.05, 1.0),
	# Covariance parameters associated with the Gaussian process
	ρ  = Uniform(2.0, 10.0),
	ν  = Uniform(0.5,  3.0),
	# Parameters of the Subbotin distribution
	μ  = Uniform(-0.5, 0.5),
	τ  = Uniform(0.3, 0.9),
	δ₁ = Uniform(1.3, 3.0)
)
parameter_names = String.(collect(keys(Ω)))

N = 16
S = expandgrid(1:N, 1:N) .- [0.5 0.5]
S = Float64.(S)
D = pairwise(Euclidean(), S, S, dims = 1)
s₀ = isodd(N) ? N/2 : (N-1)/2
s₀ = [s₀, s₀]'
s₀ = Float64.(s₀)
u = standardLaplacequantile(0.975)
h = map(norm, eachslice(S .- s₀, dims = 1))
h = Float64.(h)

# Compute versions of the above vectors/matrices that exclude s₀
s₀_idx = findfirst(x -> x == 0.0, h)
S̃  = S[1:end .!= s₀_idx, :]
h̃  = h[1:end .!= s₀_idx]

ξ = (
	Ω = Ω, p = length(Ω),
	S = S, D = D, h = h, s₀ = s₀, u = u,
	parameter_names = parameter_names,
	ρ_idx = findfirst(parameter_names .== "ρ"),
	ν_idx = findfirst(parameter_names .== "ν"),
	invtransform = x -> x^3,
	s₀_idx = s₀_idx, S̃ = S̃, h̃ = h̃
 )
