using NeuralEstimators
using SpatialDeepSets
using DataFrames
using CSV
include(joinpath(pwd(), "src/RedSea/Bootstrap.jl"))

for data_type ∈ ("regular", "irregular")

    data_path = joinpath(pwd(), "data/RedSea/" * data_type * "/")
	blocks = loadwithoutdict(data_path * "/blocks.rda", "blocks")
	region = loadwithoutdict(data_path * "/region_id.rda", "region_id")
	u      = loadwithoutdict(data_path * "/u.rda", "u")
    Z      = loadwithoutdict(data_path * "extreme_data_subset_LaplaceScale.rda", "extreme_data_L")

	# The number of fields in each bootstrap data set is variable (some groups
	# are larger than others), so we cannot store the data as a 3D array. We
	# have to store it as a vector of 2D arrays, where the number of columns varies.
    Z̃ = vcat([nonparametricbootstrap(identity, Z, blocks, B = 1, use_gpu = false) for b ∈ 1:B₂]...)

	# Compute the proprtion of exceedances in each region
	π̃ = bootstrap_threshold_exceedances(Z̃, region, u)

    CSV.write(joinpath("intermediates/RedSea/" * data_type, "empirical_exceedances_bootstrap.csv"), π̃)
end
