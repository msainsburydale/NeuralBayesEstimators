using ArgParse
arg_table = ArgParseSettings()
@add_arg_table arg_table begin
	"--arch"
		  help = "The architecture; 'CNN' or 'DNN'"
		  arg_type = String
		  required = true
      "--data_type"
            help = "The kind of data: regular or irregular"
            arg_type = String
            default = "regular"
	  "--quick"
			help = "A flag controlling whether or not a computationally inexpensive run should be done."
			action = :store_true
end
parsed_args = parse_args(arg_table)
data_type = parsed_args["data_type"]
quick     = parsed_args["quick"]
arch      = parsed_args["arch"]
@assert arch ∈ ("CNN", "DNN")
@assert !(data_type == "irregular" && arch == "CNN") "CNN cannot be used with irregular data"
@info "Red Sea study: $arch architecture for $data_type data"

intermediates_path = "intermediates/RedSea/$arch"
if arch == "DNN" intermediates_path = "intermediates/RedSea/DNN" * data_type end
if !isdir(intermediates_path) mkpath(intermediates_path) end

using NeuralEstimators
using NeuralEstimators: subsetdata
using SpatialDeepSets
using CSV
using DataFrames
using Tables
using Random: seed!
include(joinpath(pwd(), "src/RedSea/Parameters.jl"))
include(joinpath(pwd(), "src/RedSea/$arch/Simulation.jl"))
include(joinpath(pwd(), "src/RedSea/$arch/Architecture.jl"))

# ---- Sample parameters and simulation training data ----

params_path = joinpath(pwd(), "intermediates/RedSea/" * data_type * "/parameter_configurations/")
θ_train = Parameters(params_path * "train_", data_type, J = 10)
θ_val   = Parameters(params_path * "val_", data_type, J = 10)

seed!(1)
if arch == "DNN"
	ψ, ϕ = architecture(ξ.n, ξ.p)
elseif arch == "CNN"
	ψ, ϕ = architecture(ξ.p)
end
θ̂ = DeepSet(ψ, ϕ, a = "mean")


# ---- Simulate data ----

m = [1, 10, 30, 75, 150]

@info "Simulating validation data..."
sim_time = @elapsed Z_val = simulate(θ_val, ξ, maximum(m))
@info "Simulating training data..."
sim_time += @elapsed Z_train = simulate(θ_train, ξ, 2 * maximum(m))
CSV.write("$intermediates_path/sim_time.csv", Tables.table([sim_time]), header = false)


# ---- Train the neural estimator ----

batchsize = [128, 128, 64, 32, 16]

for i ∈ eachindex(m)

	mᵢ = m[i]

	@show mᵢ

	global θ̂ = train(
		  θ̂, θ_train, θ_val, Z_train, subsetdata(Z_val, 1:mᵢ), batchsize = batchsize[i],
		  savepath = "$intermediates_path/runs_ND_" * "m$(mᵢ)",
		  stopping_epochs = mᵢ == 1 ? 10 : 4
	)
end
