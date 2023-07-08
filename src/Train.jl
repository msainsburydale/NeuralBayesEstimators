# Read command line arguments
using ArgParse
arg_table = ArgParseSettings()
@add_arg_table arg_table begin
	"--model"
		help = "A relative path to the folder of the assumed model; this folder should contain scripts for defining the parameter configurations in Julia and for data simulation."
		arg_type = String
		required = true
	"--quick"
		help = "A flag controlling whether or not a computationally inexpensive run should be done."
		action = :store_true
	"--deep"
		help = "A flag controlling whether to aggregate information within the exchangeable network rather than at the end."
		action = :store_true
	"--m"
		# help = "The number of field replicates per parameter configuration during training. Providing a negative integer means a variable number of fields are used during training (i.e., --m=-30 means that between 1 and 30 fields are used per parameter configuration, and these are randomly selected and change each epoch)"
		help = "The sample size to use during training. If multiple samples sizes are given as a vector, multiple neural estimators will be trained."
		arg_type = String
end
parsed_args = parse_args(arg_table)
model = parsed_args["model"]
quick = parsed_args["quick"]
deep  = parsed_args["deep"]

# See https://stackoverflow.com/a/61297773
# m     = parsed_args["m"]; if m < 0 m = 1:(-m) end
m = let expr = Meta.parse(parsed_args["m"])
    @assert expr.head == :vect
    Int.(expr.args)
end

@show model
@show m

using NeuralEstimators
using NeuralEstimators: subsetdata
using SpatialDeepSets
using CSV
using DataFrames
using Tables
using Random: seed!

include(joinpath(pwd(), "src/$model/Parameters.jl"))
include(joinpath(pwd(), "src/$model/Simulation.jl"))
include(joinpath(pwd(), "src/Architecture.jl"))

if !isdir("intermediates/$model") mkpath("intermediates/$model") end


batchsize = 32

if model == "GaussianProcess/nuFixed"

	params_path = joinpath(pwd(), "intermediates/GaussianProcess/nuFixed/parameter_configurations/")
	θ_train     = Parameters(params_path * "train_")
	θ_val       = Parameters(params_path * "val_")

else

	K_train = 10_000
	K_val   = K_train ÷ 5

	if quick
		K_train = K_train ÷ 100
		K_val   = K_val ÷ 100
	end

	@info "Sampling validation parameters..."
	sample_time = @elapsed θ_val = Parameters(K_val, ξ)
	@info "Sampling training parameters..."
	sample_time += @elapsed θ_train = Parameters(K_train, ξ)
	CSV.write("intermediates/$model/sample_time.csv", Tables.table([sample_time]), header = false)
end

seed!(1)
ψ, ϕ = architecture(ξ.p, 0)

if deep
	θ̂ = DeepSet(ψ, ϕ)
else
	θ̂ = DeepSet(Chain(ψ, ϕ), identity)
end

savepath = "intermediates/$model/runs_N$(deep ? "D" : "")"

if length(m) == 1

	m = m[1]

	θ̂ = train(
		  θ̂, θ_train, θ_val, simulate, m = m, batchsize = batchsize,
		  savepath = savepath * "_m$(m)",
		  stopping_epochs = m == 1 ? 10 : 4
	)
else

	@info "Simulating validation data..."
	sim_time = @elapsed Z_val = simulate(θ_val, maximum(m))
	@info "Simulating training data..."
	sim_time += @elapsed Z_train = simulate(θ_train, 2 * maximum(m))
	CSV.write("intermediates/$model/sim_time.csv", Tables.table([sim_time]), header = false)

	for mᵢ ∈ m

		@show mᵢ

		optimiser = mᵢ == 1 ? ADAM(1e-4) : ADAM(1e-4/10)

		global θ̂ = train(
			  θ̂, θ_train, θ_val, Z_train, subsetdata(Z_val, 1:mᵢ), batchsize = batchsize,
			  savepath = savepath * "_m$(mᵢ)", optimiser = optimiser,
			  stopping_epochs = mᵢ == 1 ? 10 : 4
		)
	end

end
