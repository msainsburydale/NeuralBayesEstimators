using ArgParse
arg_table = ArgParseSettings()
@add_arg_table arg_table begin
	"--arch"
		  help = "The architecture; 'CNN' or 'FC'"
		  arg_type = String
		  required = true
      "--data_type"
            help = "The kind of data: regular or irregular"
            arg_type = String
            default = "regular"
end
parsed_args = parse_args(arg_table)
data_type = parsed_args["data_type"]
arch      = parsed_args["arch"]
@assert arch ∈ ("CNN", "FC")
@assert !(data_type == "irregular" && arch == "CNN") "CNN cannot be used with irregular data"

intermediates_path = "intermediates/RedSea/$arch"
if arch == "FC" intermediates_path = "intermediates/RedSea/FC" * data_type end

using SpatialDeepSets
using Distances: pairwise, Euclidean
using LinearAlgebra
using DataFrames
using CSV
using Random: seed!
include(joinpath(pwd(), "src/RedSea/Parameters.jl"))
include(joinpath(pwd(), "src/RedSea/$arch/Simulation.jl"))
include(joinpath(pwd(), "src/RedSea/$arch/Architecture.jl"))
params_path = joinpath(pwd(), "intermediates/RedSea/" * data_type * "/parameter_configurations/")
θ_test      = Parameters(params_path * "test_", data_type)
θ_scenarios = Parameters(params_path * "scenarios_", data_type)

relative_loadpath = intermediates_path
relative_savepath = intermediates_path * "/Estimates"
savepath = joinpath(pwd(), relative_savepath)
if !isdir(savepath) mkdir(savepath) end


# ---- Load neural estimators ----

# Don't need to worry about constructing a piecewise estimator because we are
# estimating with fixed m ≈ 140.

if arch == "CNN"
	NNs = loadneuralestimators(joinpath(pwd(), relative_loadpath), architecture, ξ.p)
elseif arch == "FC"
	NNs = loadneuralestimators(joinpath(pwd(), relative_loadpath), architecture, ξ.p, ξ.n)
end
networks = NNs.estimators
titles = NNs.titles

# Select the desired estimator.
ND_titles = titles[occursin.("ND", titles)]
m = maximum(parse.(Int64, replace.(ND_titles, "ND_m" => "")))
title = "ND_m$m"
net = networks[titles .== title][1]
net = net |> gpu # move to the GPU


# ---- Estimates from real data: Fitted values ----

parameter_names = ξ.parameter_names

# number of replicates available at the estimation stage, mₑ:
mₑ = loadwithoutdict(joinpath(pwd(), "data/RedSea/" * data_type *"/m_e.rda"), "m_e") |> Int

# Load the data
data_path = joinpath(pwd(), "data/RedSea/" * data_type * "/")
blocks   = loadwithoutdict(data_path * "/blocks.rda", "blocks")
region   = loadwithoutdict(data_path * "/region_id.rda", "region_id")
u        = loadwithoutdict(data_path * "/u.rda", "u")
Z_RedSea = loadwithoutdict(data_path * "extreme_data_subset_LaplaceScale.rda", "extreme_data_L")
Z_RedSea = Float32.(Z_RedSea)
Z_RedSea = cbrt.(Z_RedSea) # variance stablising transformation

if arch == "FC"
	Z_RedSea = reshape(Z_RedSea, size(Z_RedSea, 1), 1, size(Z_RedSea, 2))
elseif arch == "CNN"
	# Pad the grid with zeros, so that it can be processed by the CNN
	# Information needed to pad the grid with zeros.
	# This step also converts the data from an n x mₑ matrix to a hᵢ x wᵢ x 1 x mₑ
	# array, where hᵢ and wᵢ are the height and width, respectively, of the full
	# grid (i.e., the grid with zeros including over the padding region).
	Z_RedSea_matrix = map(1:mₑ) do i
		A = fill(0f0, ξ.height, ξ.width)
		A[ξ.data_idx] = Z_RedSea[:, i]
		A
	end
	Z_RedSea = cat(Z_RedSea_matrix..., dims = 4)
end

# move to the GPU
Z_RedSea = Z_RedSea|> gpu

# Fit the parameters
@elapsed net([Z_RedSea]) # one run to compile
t = @elapsed θ̂ = net([Z_RedSea])
θ̂ = θ̂ |> cpu # move back to the CPU
CSV.write(joinpath(savepath, "real_data_est_time.csv"), Tables.table([t]), header=false)
CSV.write(joinpath(savepath, "real_data_estimates.csv"), DataFrame(θ̂', parameter_names))

# Create a Parameters object from the fitted parameters
fitted_params = Parameters(θ̂, ξ, J = 1)

# Simulate from the fitted model for plotting
Z_df = simulatedataframe(fitted_params, ξ)
CSV.write(joinpath(savepath, "fitted_model_simulations.csv"), Z_df)

# ---- Bootstrapping ----

include(joinpath(pwd(), "src/RedSea/Bootstrap.jl"))

# Bootstrap samples of θ̂ (non-parametric)
blocks = loadwithoutdict(joinpath(pwd(), "data/RedSea/" * data_type * "/blocks.rda"), "blocks")
nonparametricbootstrap(net, Z_RedSea, blocks, B = 10) # one run to compile the function
t = @elapsed θ̃ = nonparametricbootstrap(net, Z_RedSea, blocks, B = B)
θ̃ = θ̃ |> cpu
CSV.write(joinpath(savepath, "bootstrap_samples_nonparametric.csv"), DataFrame(θ̃', parameter_names))
CSV.write(joinpath(savepath, "bootstrap_time_nonparametric.csv"), Tables.table([t]), header=false)

# Simulate from the model using the bootstrap samples, θ̃, and compute the
# proprtion of threshold exceedances. First, construct a Parameters object for θ̃.
θ̃_params = Parameters(Float64.(θ̃), ξ, J = 1)

Z̃ = simulate(θ̃_params, ξ, B₂)
Z̃ = broadcast.(ξ.invtransform, Z̃)
if arch == "FC"
	Z̃ = dropdims.(Z̃, dims = 2)
	π̃ = bootstrap_threshold_exceedances(Z̃, region, u)
elseif arch == "CNN"
	π̃ = bootstrap_threshold_exceedances(Z̃, region, u, data_idx)
end
CSV.write(joinpath(savepath, "threshold_bootstrap_samples_nonparametric.csv"), π̃)
# threshold_exceedances(π̃)

# ---- Assessing the estimator ----

estimate(
	[net], ξ, θ_test;  m = [1, 5, 10, 20, 30, 60, 90, 120, 150],
	estimator_names = [title], parameter_names = ξ.parameter_names,
	save = [savepath, "test"]
	)

# Add θ̂ to test if the estimator performs well when θ = θ̂.
scenarios_and_fitted_params = Parameters(θ_scenarios, fitted_params)

estimate(
	[net], ξ, θ_scenarios; m = [mₑ], num_rep = 100,
	estimator_names = [title], parameter_names = ξ.parameter_names,
	save = [savepath, "scenarios"]
	)

seed!(1)
fields_df = simulatedataframe(scenarios_and_fitted_params, ξ)
CSV.write(joinpath(savepath, "fields_scenarios.csv"), fields_df)
