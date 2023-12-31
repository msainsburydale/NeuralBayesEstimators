# Read command line arguments
using ArgParse
arg_table = ArgParseSettings()
@add_arg_table arg_table begin
	"--quick"
		help = "A flag controlling whether or not a computationally inexpensive run should be done."
		action = :store_true
	"--model"
		help = "A relative path to the folder of the assumed model."
		arg_type = String
		required = true
	"--expert"
		help = "A flag controlling whether or not to use expert summary statistics in the second stage."
		action = :store_true
	"--m"
		help = "The number of field replicates per parameter configuration during training. Providing a negative integer means a variable number of fields are used during training (i.e., --m=-30 means that between 1 and 30 fields are used per parameter configuration, and these are randomly selected and change each epoch)"
		arg_type = Int
		default = -150
	"--aggregation"
		help = "The aggregation to use; can be 'sum' (the default), 'mean', or 'logsumexp'."
		arg_type = String
		default = "mean"
end
parsed_args = parse_args(arg_table)
quick  = parsed_args["quick"]
model  = parsed_args["model"]
expert = parsed_args["expert"]
aggregation = parsed_args["aggregation"]
m = parsed_args["m"]; if m < 0 m = 1:(-m) end

using NeuralEstimators
using SpatialDeepSets
using Flux: logsumexp
using Statistics: mean, sum
using CSV
using DataFrames
using Tables
using Random: seed!

include(joinpath(pwd(), "src", "Univariate", "Architecture.jl"))
include(joinpath(pwd(), "src", "Univariate", model, "Parameters.jl"))
include(joinpath(pwd(), "src", "Univariate", model, "Simulation.jl"))

m_title = typeof(m) <: AbstractRange ? "$(m.start)to$(m.stop)" : "$m"
title =  "m$(m_title)$(expert ? "_expert" : "")"
savepath = joinpath("intermediates", "Univariate", model, "runs_NN_$title")

if !isdir(savepath) mkpath(savepath) end

qₜ = 100 # dimension of the learned summary statistic
CSV.write(
	joinpath(pwd(), savepath, "qt.csv"),
	Tables.table([qₜ]), header = false
	)
CSV.write(
	joinpath(pwd(), savepath, "aggregation.csv"),
	Tables.table([aggregation]), header = false
	)


batchsize = 128
epochs = quick ? 5 : 100

# size of the the data and the number of parameters.
n = 1
p = 1

seed!(1)
if expert
	@info "Using the DeepSets representation with expert summary statistics"
	# Note that we can use maximum rather than logsumexp because we don't need
	# to differentiate through S(⋅), since S(⋅) doesn't contain any parameters.
	# S = [samplesize, maximum]
	S = [inversesamplesize]
	qₛ = length(S)
	@info "qₛ = $(qₛ)"
	ψ, ϕ = architecture(n, qₜ, qₛ, p)
	θ̂ = DeepSetExpert(ψ, ϕ, S, a = aggregation)
else
	@info "Using the DeepSets representation without expert summary statistics"
	qₛ = 0
	ψ, ϕ = architecture(n, qₜ, 0, p)
	θ̂ = DeepSet(ψ, ϕ, a = aggregation)
end

K = quick ? 1000 : 100000

θ̂ = train(
	  θ̂, Parameters, simulate, K = K, m = m, epochs = epochs, ξ = ξ,
	  savepath = savepath, batchsize = batchsize,
	  use_gpu = qₛ < 2 # getting errors when training with multiple expert statistics: fix later
)
