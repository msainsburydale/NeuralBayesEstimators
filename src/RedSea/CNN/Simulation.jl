include(joinpath(pwd(), "src", "ConditionalExtremes", "Simulation.jl"))

# Since we're using a CNN, we need to convert each field to a rectangular
# array, padding with zeros to fill out the rectangle.
function reshapeZ(Z::A, ξ)  where {T <: Number, A <: AbstractArray{T, 2}}
	return pad_grid(Z, ξ.data_idx, ξ.height, ξ.width)
end


function pad_grid(Z::A, data_idx, height, width) where {T <: Number, A <: AbstractArray{T, 2}}
	Z = map(z -> pad_grid(z, data_idx, height, width), eachcol(Z))
	Z = cat.(Z; dims = 4) # add singleton third and fourth dimensions
	Z = stackarrays(Z)    # after stacking, the third dimension remains singleton to indicate univariate data in Flux language
	return Z
end

function pad_grid(Z::A, data_idx, height, width) where {T <: Number, A <: AbstractArray{T, 1}}
	full_grid = fill(zero(T), height, width)
	full_grid[data_idx] = Z
	return full_grid
end

import SpatialDeepSets: simulatedataframe
function simulatedataframe(params::Parameters, ξ)

	S = ξ.S
	n = size(S, 1)

	θ = params.θ
	p = size(θ, 1)
	K = size(θ, 2)

	n_fields = max(binomial(p, 2), 16)

	T = eltype(θ)
	fields_df = DataFrame(Z = T[], replicate = Int[], scenario = Int[], x = T[], y = T[])
	scenario  = repeat(1:K, inner = n)

	# The coordinates of a single field
	x = S[:, 1]
	y = S[:, 2]
	# Repeat x and y for the number of fields we simulate in each loop
	x = repeat(x, outer = K)
	y = repeat(y, outer = K)

	for i ∈ 1:n_fields
		Z  = simulate(params)
		Z  = vcat([vec(Z[:, :, 1, j][ξ.data_idx]) for j ∈ 1:K]...)
		df = DataFrame(Z = Z, scenario = scenario, replicate = i, x = x, y = y)
		append!(fields_df, df)
	end

	return fields_df
end
