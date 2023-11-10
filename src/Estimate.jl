using ArgParse
arg_table = ArgParseSettings()
@add_arg_table arg_table begin
	"--model"
		help = "A relative path to the folder of the assumed model."
		arg_type = String
		required = true
	"--ML"
        help = "Flag indicating whether we should conduct likelihood estimation."
    	action = :store_true
	"--quick"
		help = "A flag controlling whether or not a computationally inexpensive run should be done."
		action = :store_true
	"--linesearch"
		help = "A flag controlling whether or not to perform the line-search for the optimal fixed cut-off distance, d, when choosing pairs for the pairwise likelihood."
		action = :store_true
end
parsed_args = parse_args(arg_table)
model = parsed_args["model"]
ML    = parsed_args["ML"]
quick = parsed_args["quick"]
linesearch = parsed_args["linesearch"]

using NeuralEstimators
using SpatialDeepSets
using DataFrames
using CSV
using Random: seed!

relative_loadpath = joinpath("intermediates", model)
relative_savepath = joinpath("intermediates", model, "Estimates")
savepath = joinpath(pwd(), relative_savepath)
if !isdir(savepath) mkdir(savepath) end

include(joinpath(pwd(), joinpath("src", model, "Parameters.jl")))
include(joinpath(pwd(), joinpath("src", model, "Simulation.jl")))
include(joinpath(pwd(), joinpath("src", "Architecture.jl")))

if model == "GaussianProcess/nuFixed"
	params_path = joinpath(pwd(), "intermediates", model, "parameter_configurations")
	θ_test      = Parameters(joinpath(params_path, "test_"))
	θ_scenarios = Parameters(joinpath(params_path, "scenarios_"))
else
	seed!(1)
	K = quick ? 10 : 500
	θ_test = Parameters(K, ξ)

	# Sample randomly from Ω to generate a small set of parameters used for
	# visualising the joint distributions of the estimators.
	seed!(1)
	θ_scenarios = Parameters(5, ξ, J = 1)

end

println("\nStarting Estimation for $model model...")
println("Number of CPU threads available for data simulation: $(Threads.nthreads())")

# ---- Load neural estimators ----

loadpath   = joinpath(pwd(), relative_loadpath)
NNs        = loadneuralestimators(loadpath, architecture, ξ.p)
estimators = NNs.estimators
titles     = NNs.titles

# ---- Piecewise neural estimator ----

# Construct a piecewise neural estimator from the Deep Set neural estimators
# trained for specific values of m (the sample size). Note that m needs to match
# the values used during training.
m = [1,10,30,75,150]
m_cutoffs = [1, 20, 50, 100]
@assert length(m_cutoffs) == length(m) - 1
target_order = "ND_m" .* string.(m)
indices      = indexin(target_order, titles)
θ̂_piecewise  = PiecewiseEstimator(estimators[indices], m_cutoffs)

estimators = [estimators..., θ̂_piecewise]
push!(titles, "ND_piecewise")


# ---- Select the neural estimators we wish to estimate with ----

neural_estimators = ["N_m1", "ND_piecewise"]
idx = [findfirst(occursin.(estimator, titles)) for estimator ∈ neural_estimators]
titles = titles[idx]
estimators = estimators[idx]

println("Neural network estimators: $(join(titles, ", "))")


# ---- Likelihood functions ----

if ML
	println("Number of CPU threads available for likelihood estimation: $(Threads.nthreads())")
	include(joinpath(pwd(), "src/$model/ML.jl"))

	if model == "Schlather"
		PL3(y, ξ) = PL(y, ξ, 3)
		likelihood_estimators = [PL3]
		likelihood_titles = ["PL3"]
	else
		likelihood_estimators = [likelihoodestimator]
		likelihood_titles = ["ML"]
	end

	println("Likelihood estimators: $(join(likelihood_titles, ", "))")
end


# ---- Estimation over the test set ----

all_m   = [1, 5, 10, 20, 30, 60, 90, 120, 150]
estimates = estimate(
	estimators, ξ, θ_test; m = all_m,
	estimator_names = titles, parameter_names = ξ.parameter_names,
	save = [savepath, "test"]
)

if ML

	# NB The following code is needed to generate Figure S12 of the
	#    Supplementary Material (which illustrates the line search for the fixed
	#    cut-off distance, d). We do not run it by default (it takes a long time
	#    to run and is not central to the paper), but this can be changed by
	#    adding the flag --linesearch in the controlling .sh file.

	if model == "Schlather" && linesearch
		# Estimate using several different pairwise likelihood estimators to
		# show that d = 3 is optimal:
		K = quick ? 5 : 100
		θ_testsmall = Parameters(K, ξ)

		PL2(y, ξ) = PL(y, ξ, 2)
		PL3(y, ξ) = PL(y, ξ, 3)
		PL6(y, ξ) = PL(y, ξ, 6)
		PL9(y, ξ) = PL(y, ξ, 9)
		PLF(y, ξ) = PL(y, ξ, Inf)

		estimates = estimate(
			[PL2, PL3, PL6, PL9, PLF], (ξ..., θ₀ = θ_testsmall.θ), θ_testsmall;
			m = [1, 10, 20, 30, 60, 100, 150], use_ξ = true, use_gpu = false,
			estimator_names = ["PL2", "PL3", "PL6", "PL9", "PLF"],
			parameter_names = ξ.parameter_names,
			save = [savepath, "likelihood_test_allpairwise"]
			)
	end

	ξ₀ = (ξ..., θ₀ = θ_test.θ)

	estimates = estimate(
		likelihood_estimators, ξ₀, θ_test;
		m = all_m, use_ξ = true, use_gpu = false,
		estimator_names = likelihood_titles,
		parameter_names = ξ.parameter_names,
		save = [savepath, "likelihood_test"]
		)
end


# ---- Estimation over several "scenarios", a subset of the test set ----

all_m = [1, 30, 150]

num_rep = 100

estimates = estimate(
	estimators, ξ, θ_scenarios;
	m = all_m,
	num_rep = num_rep,
	estimator_names = titles,
	parameter_names = ξ.parameter_names,
	save = [savepath, "scenarios"]
)

if ML

	ξ₀ = (ξ..., θ₀ = θ_scenarios.θ)

	estimates = estimate(
		likelihood_estimators, ξ₀, θ_scenarios; m = all_m,
		num_rep = quick ? 5 : num_rep, # likelihood estimation is slow, so only do a few if quick=true
		use_ξ = true, use_gpu = false,
		estimator_names = likelihood_titles,
		parameter_names = ξ.parameter_names,
		save = [savepath, "likelihood_scenarios"]
	)

end
