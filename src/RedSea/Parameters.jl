using SpatialDeepSets
import SpatialDeepSets: Parameters
using SpatialDeepSets: objectindices
using LinearAlgebra
using Random: seed!

# The only real difference between here and the simulation study version is that
# we load the matrix S (since it corresponds to a physical domain now) and the
# extremal threshold, u, which is a quantity defined in terms of the data.
# There's also no prior measure defined here, since we sampled the parameters in
# a separate R script to reduce the computational cost (they are needed in
# multiple studies).

data_path = "data/RedSea/" * data_type

S      = loadwithoutdict(data_path * "/S.rda", "S")
D      = loadwithoutdict(data_path * "/D.rda", "D")
s₀     = loadwithoutdict(data_path * "/s0.rda", "s0")
h      = map(norm, eachslice(S .- s₀, dims = 1))
s₀_idx = findfirst(x -> x == 0.0, h)
u      = loadwithoutdict(data_path * "/u.rda", "u")

# Compute versions of the above vectors/matrices that exclude s₀
s₀_idx = findfirst(x -> x == 0.0, h)
S̃  = S[1:end .!= s₀_idx, :]
h̃  = h[1:end .!= s₀_idx]

parameter_names = ["κ", "λ", "β", "ρ", "ν", "μ", "τ", "δ₁"]

ξ = (
	p = length(parameter_names),
	n = size(S, 1),
	S = S, D = D, h = h, s₀ = s₀, u = u,
	parameter_names = parameter_names,
	ρ_idx = findfirst(parameter_names .== "ρ"),
	ν_idx = findfirst(parameter_names .== "ν"),
	invtransform = x -> x^3,
	s₀_idx = s₀_idx, S̃ = S̃, h̃ = h̃
 )


 if data_type == "regular"

	 # Between the regular and irregular Red Sea applications, the only difference
	 # here is that, for the regular data, we need to load the height and width of
	 # the padded data array. We also need data_idx, where data_idx[i] gives the
	 # linear index with respect to the full grid corresponding to the spatial
	 # location, S[i, ]. These objects are stored in the invariant model
	 # information, ξ.

	 width     = loadwithoutdict(data_path * "/width.rda", "width")
	 height    = loadwithoutdict(data_path * "/height.rda", "height")
	 data_idx  = loadwithoutdict(data_path * "/data_idx.rda", "data_idx")

	 ξ = (ξ..., width = width, height = height, data_idx = data_idx)
 end


function Parameters(params_path::String, data_type::String; J::Integer = 1)

 	data_path = "data/RedSea/" * data_type

	# transpose() because we want to convert from long to wide format (Flux
	# assumes that the last dimension of an array is the observation dimension),
	# and copy() forces a regular array rather than an object of type Transpose.
	θ     = loadwithoutdict(params_path * "xi.rda", "xi") |> transpose |> copy
	chols = loadwithoutdict(params_path * "chols.rda", "chols")

	# Convert θ to Float32 for more efficient training.
	θ = Float32.(θ)
	θ = repeat(θ, inner = (1, J))

	return Parameters(θ, chols, objectindices(chols, θ))
end



# function Parameters(θ̂::Matrix, params::Parameters, ξ)
#
# 	ρ = θ̂[ξ.ρ_idx, :]
# 	ν = θ̂[ξ.ν_idx, :]
# 	D = Float64.(ξ.D)
# 	L = Folds.map(ρ, ν) do ρᵢ, νᵢ LinearAlgebra.cholesky(matern.(D, ρᵢ, νᵢ)).L end
# 	L = cat(L..., dims = 3) # NB: splatting within cat is slow if there are a large number of terms to splat; if this becomes an issue, just use stackarrays()
# 	L = Float32.(L)
# 	chols = reshape(L, (size(L)..., 1))
#
# 	return Parameters(θ = θ̂, chols = chols, S = params.S, D = params.D, s₀ = params.s₀, u = params.u)
# end
# function Parameters(θ̂, params::Parameters)
#
# 	param_names = params.param_names
# 	ρ_idx = findall(contains.(param_names, "ρ"))[1]
# 	ν_idx = findall(contains.(param_names, "ν"))[1]
# 	ρ = exp.(θ̂[ρ_idx, :])
# 	ν = exp.(θ̂[ν_idx, :])
# 	D = Float64.(params.D)
# 	L = Folds.map(ρ, ν) do ρᵢ, νᵢ LinearAlgebra.cholesky(matern.(D, ρᵢ, νᵢ)).L end
# 	L = cat(L..., dims = 3) # NB: splatting within cat is slow if there are a large number of terms to splat; if this becomes an issue, just use the stack() function I use elsewhere.
# 	L = Float32.(L)
# 	chols = reshape(L, (size(L)..., 1))
#
# 	return Parameters(θ = θ̂, chols = chols, S = params.S, D = params.D, s₀ = params.s₀, u = params.u)
# end
#
# # ---- constructor used by _getparameters() ----
#
# function Parameters(params::Parameters; θ, chols)
# 	Parameters(
# 	θ = θ,
#     chols = chols,
# 	S  = params.S,
# 	D  = params.D,
# 	s₀ = params.s₀,
# 	u  = params.u
# 	)
# end
