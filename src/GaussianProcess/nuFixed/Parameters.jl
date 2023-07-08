using SpatialDeepSets
import SpatialDeepSets: Parameters
using Distances: pairwise, Euclidean
using LinearAlgebra
using Distributions: Uniform
using Random: seed!

# NB The parameters generated in R need to be contained within these bounds
Ω = (
	σ = Uniform(0.001, 1.5),
	ρ = Uniform(1.0, 30.0)
)
parameter_names = String.(collect(keys(Ω)))
S = expandgrid(1:16, 1:16)
S = Float64.(S)
D = pairwise(Euclidean(), S, S, dims = 1)
ξ = (
	Ω = Ω, S = S, D = D, p = length(Ω),
	parameter_names = parameter_names,
	ρ_idx = findfirst(parameter_names .== "ρ"),
	invtransform = identity
)

function Parameters(path::String)

	# transpose() because we want to convert from long to wide format (Flux
	# assumes that the last dimension of an array is the observation dimension),
	# and copy() forces a regular array rather than an object of type Transpose.
	θ        = loadwithoutdict(path * "xi.rda", "xi") |> transpose |> copy
	S        = loadwithoutdict(dirname(path) * "/S.rda", "S")
	D        = loadwithoutdict(dirname(path) * "/D.rda", "D")
	chols    = loadwithoutdict(path * "chols.rda", "chols")

	K = size(θ, 2)
	num_chols = size(chols, 3)
	num_reuse = K ÷ num_chols
	chol_idx  = repeat(1:num_chols, inner = num_reuse)

	νFixed = length(unique(θ[3, :])) == 1
	if νFixed θ = θ[1:2, :] end

	return Parameters(Float32.(θ), Float64.(chols), chol_idx)
end

function Parameters(K::Integer, ξ; J::Integer = 10)

	# All parameters not associated with the Gaussian process
	σ = rand(ξ.Ω.σ, K * J)

	# Covariance parameters associated with the Gaussian process
	ρ = rand(ξ.Ω.ρ, K)
	ν = fill(1.0, K)
	chols = maternchols(ξ.D, ρ, ν)
	ρ = repeat(ρ, inner = J)

	# Concatenate into a matrix and convert to Float32 for efficiency
	@assert ξ.ρ_idx == 2 "The parameter ρ should be stored in the second row of θ"
	θ = hcat(σ, ρ)'
	θ = Float32.(θ)

	Parameters(θ, chols, objectindices(chols, θ))
end
