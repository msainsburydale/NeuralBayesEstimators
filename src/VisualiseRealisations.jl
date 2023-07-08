# Read command line arguments
using ArgParse
arg_table = ArgParseSettings()
@add_arg_table arg_table begin
	"--model"
		help = "A relative path to the folder of the assumed model."
		arg_type = String
		required = true
end
parsed_args = parse_args(arg_table)
model = parsed_args["model"]

println("\nVisualising realisations from the $model model...")

using NeuralEstimators
using SpatialDeepSets
using DataFrames
using CSV
using Random: seed!

include(joinpath(pwd(), "src/$model/Parameters.jl"))
include(joinpath(pwd(), "src/$model/Simulation.jl"))
if model == "GaussianProcess/nuFixed"
	params_path = joinpath(pwd(), "intermediates/GaussianProcess/nuFixed/parameter_configurations/")
	θ_test      = Parameters(params_path * "test_")
	θ_scenarios = Parameters(params_path * "scenarios_")
else
	# Sample randomly from Ω to generate a small set of parameters used for
	# visualising the joint distributions of the estimators.
	seed!(1)
	θ_scenarios = Parameters(5, ξ, J = 1)
end

seed!(1)
fields_df = simulatedataframe(θ_scenarios, ξ)
fields_df[:, :Z] = ξ.invtransform.(fields_df[:, :Z])
Z = fields_df.Z
println("Field simulations: Minimum = $(minimum(Z)), Mean = $(sum(Z)/length(Z)), Maximum = $(maximum(Z))\n")

# Save the fields
savepath = joinpath(pwd(), "intermediates/$model")
!ispath(savepath) && mkpath(savepath)
CSV.write(joinpath(savepath, "fields.csv"), fields_df)

# Also save the parameters
θ_df = DataFrame(θ_scenarios.θ', ξ.parameter_names)
CSV.write(joinpath(savepath, "scenario_parameters.csv"), θ_df)
