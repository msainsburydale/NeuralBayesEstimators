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

using NeuralEstimators
using SpatialDeepSets
using DataFrames
using CSV
using Random: seed!

relative_loadpath = joinpath("intermediates", "Univariate", model)
savepath = joinpath(pwd(), relative_loadpath, "Estimates")
if !isdir(savepath) mkdir(savepath) end

include(joinpath(pwd(), "src", "Univariate", model, "BayesEstimator.jl"))
include(joinpath(pwd(), "src", "Univariate", model, "Parameters.jl"))
include(joinpath(pwd(), "src", "Univariate", model, "Simulation.jl"))
include(joinpath(pwd(), "src", "Univariate", "Architecture.jl"))

seed!(1)
θ_test = Parameters(1000, ξ)

# ---- Load NN estimators ----

# (Hard coded for now) size of the the data and the number of parameters.
n = 1
p = 1

"""
	loadpath: String giving the path containing folders called 'runs_*',
	where * denotes the name of a NN (e.g., runs_N30).
"""
function loadneuralestimators(loadpath)

	# Network names
	titles = readdir(loadpath)
	titles = titles[findall(occursin.("runs_", titles))]
	titles = replace.(titles, "runs_" => "")
	@info "Neural estimators: $(join(titles, ", "))"

	estimators = map(titles) do title

		# Get the size of the latent space and the aggregation used
		path = joinpath(pwd(), loadpath, "runs_$title")
		qₜ = CSV.read(joinpath(path, "qt.csv"), DataFrame, header = false)[1, 1]
		aggregation = CSV.read(joinpath(path, "aggregation.csv"), DataFrame, header = false)[1, 1]
		aggregation = String(aggregation)

		# Load the network
		if contains(title, "expert")
			S = [inversesamplesize]
			qₛ = length(S)
			ψ, ϕ = architecture(n, qₜ, qₛ, p)
			network = DeepSetExpert(ψ, ϕ, S, a = aggregation)
		else
			qₛ = 0
			ψ, ϕ = architecture(n, qₜ, qₛ, p)
			network = DeepSet(ψ, ϕ, a = aggregation)
		end

		Flux.loadparams!(network, loadbestweights(path))

		return network
	end

	return (estimators = estimators, titles = titles)
end


loadpath   = joinpath(pwd(), relative_loadpath)
NNs        = loadneuralestimators(loadpath)
estimators = NNs.estimators
titles     = NNs.titles
use_gpu = fill(true, length(estimators))
use_ξ   = fill(false, length(estimators))

# ---- Analytic estimators ----

if @isdefined MLE
	@info "Including the maximum likelihood estimator"
	estimators = [estimators..., MLE]
	titles = [titles..., "ML"]
	push!(use_gpu, false)
	push!(use_ξ,   true)
end

if @isdefined BayesEstimator
	@info "Including the Bayes estimator"
	estimators = [estimators..., BayesEstimator]
	titles = [titles..., "BayesEstimator"]
	push!(use_gpu, false)
	push!(use_ξ,   true)
end

if @isdefined MAPEstimator
	@info "Including the MAP estimator"
	estimators = [estimators..., MAPEstimator]
	titles = [titles..., "MAPEstimator"]
	push!(use_gpu, false)
	push!(use_ξ,   true)
end

if @isdefined OAATEstimator
	@info "Including the one-at-a-time estimator"
	estimators = [estimators..., OAATEstimator]
	titles = [titles..., "OAATEstimator"]
	push!(use_gpu, false)
	push!(use_ξ,   true)
end

# ---- Estimation ----

if model == "Uniform"

	all_m   = [10]
	num_rep = 30_000
	estimate(
	estimators, ξ, θ_scenarios, m = all_m, num_rep = num_rep,
	estimator_names = titles, use_gpu = use_gpu, save = [savepath, "scenarios"],
	parameter_names = ξ.parameter_names, use_ξ = use_ξ
	)

else
	all_m   = [1, 5, 15, 30, 60, 90, 120, 150]

	num_rep = 10
	estimate(
	estimators, ξ, θ_test, m = all_m, num_rep = num_rep, use_gpu = use_gpu,
	estimator_names = titles, save = [savepath, "test"],
	parameter_names = ξ.parameter_names, use_ξ = use_ξ
	)

	num_rep = 10_000
	estimate(
	estimators, ξ, θ_scenarios, m = all_m, num_rep = num_rep, use_gpu = use_gpu,
	estimator_names = titles, save = [savepath, "scenarios"],
	parameter_names = ξ.parameter_names, use_ξ = use_ξ
	)
end

# Test code:
# m = 10
# y = simulate(θ_scenarios, ξ, m, num_rep)
# θ  = θ_scenarios.θ
# θ̂₁ = MLE(y, ξ)
# θ̂₂ = BayesEstimator(y, ξ)
# θ̂₃ = OAATEstimator(y, ξ)
