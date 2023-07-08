using SpatialDeepSets: loadwithoutdict
using Statistics: mean

B  = 400  # number of bootstrap samples of θ̂
B₂ = 1000 # number of fields used for each bootstrap threshold exceedance estimate

function proportion_threshold_exceedances(z̃, region, u)

	@assert length(region) == size(z̃, 1) "The length of `region` (a grouping vector) does not match the number of spatial locations, `size(z̃, 1)`."

	z̃_df = DataFrame(hcat(region, z̃), :auto)
	rename!(z̃_df, :x1 => :region)
	z̃_df = stack(z̃_df, 2:ncol(z̃_df))  # equivalent of melt() in R, treating region (the first column) as the id variable
	z̃_df = groupby(z̃_df, :region)     # group the data by region
	z̃_df = combine(z̃_df, :value => (Z -> mean(Z .> u)) => :prop_exceedances) # equivalent of summarise in R

	return z̃_df
end

function bootstrap_threshold_exceedances(Z̃::V, region, u) where {T <: Number, A <: AbstractArray{T, 2}, V <: AbstractArray{A}}
	π̃ = DataFrame(region = Float64[], prop_exceedances = Float64[], b = Int[])
	for b ∈ eachindex(Z̃)
		df = proportion_threshold_exceedances(Z̃[b], region, u)
		df[!, :b] = repeat([b], nrow(df))
		π̃ = vcat(π̃, df)
	end
	return π̃
end


function bootstrap_threshold_exceedances(Z̃::V, region, u, data_idx) where {T <: Number, A <: AbstractArray{T, 4}, V <: AbstractArray{A}}

	Z̃_reshaped = map(Z̃) do z̃
		z̃ = dropdims(z̃, dims = 3)
		z̃ = mapslices(x -> vec(x[data_idx]), z̃, dims = [1, 2])
		z̃ = dropdims(z̃, dims = 2)
		z̃
	end

	return bootstrap_threshold_exceedances(Z̃_reshaped, region, u)
end

using Statistics: quantile
function threshold_exceedances(π̃::DataFrame)
	π̃ = groupby(π̃, :region)
	π̃ = combine(π̃,
				:prop_exceedances => mean => :mean,
				:prop_exceedances => (x->quantile(x, 0.025)) => :lower,
				:prop_exceedances => (x->quantile(x, 0.975)) => :upper)
	return π̃
end
