# Read command line arguments
using ArgParse
arg_table = ArgParseSettings()
@add_arg_table arg_table begin
	"--quick"
		help = "A flag controlling whether or not a computationally inexpensive run should be done."
		action = :store_true
end
parsed_args = parse_args(arg_table)
quick = parsed_args["quick"]

using NeuralEstimators
using NeuralEstimators: subsetdata
using SpatialDeepSets
using Random: seed!

model = "GaussianProcess/nuVaried"
model = replace(model, "/" => PATH_SEPARATOR)
include(joinpath(pwd(), "src", model, "Parameters.jl"))
include(joinpath(pwd(), "src", model, "Simulation.jl"))
include(joinpath(pwd(), "src", "Architecture.jl"))

K_train = 10_000
K_val   = K_train ÷ 5
if quick
	K_train = K_train ÷ 100
	K_val   = K_val ÷ 100
end
epochs = quick ? 2 : 300
stopping_epochs = 3
batchsize = 32

m = 30

@info "Sampling validation parameters..."
θ_val = Parameters(K_val, ξ)
@info "Simulating validation data..."
Z_val = simulate(θ_val, m)

@info "Sampling training parameters..."
θ_train = Parameters(K_train, ξ)
@info "Simulating training data..."
Z_train = simulate(θ_train, 5m)



# ---- Train the NN₃₀ with pre-training ----

@info "Training with pre-training"

seed!(1)
ψ, ϕ = architecture(ξ.p, 0)
θ̂ = DeepSet(ψ, ϕ)

θ̂ = train(
	  θ̂, θ_train, θ_val, Z_train, subsetdata(Z_val, 1:1),
	  batchsize = batchsize, epochs = epochs,
	  savepath = joinpath("intermediates", "Pretraining", "Pretrained", "runs_N1"),
	  stopping_epochs = stopping_epochs
)

θ̂ = train(
	  θ̂, θ_train, θ_val, Z_train, Z_val,
	  batchsize = batchsize, epochs = epochs, optimiser = ADAM(1e-5),
	  savepath = joinpath("intermediates", "Pretraining", "Pretrained", "runs_N30"),
	  stopping_epochs = stopping_epochs
)

# ---- Train the NN₃₀ without pre-training ----

@info "Training without pre-training"

seed!(1)
ψ, ϕ = architecture(ξ.p, 0)
θ̂ = DeepSet(ψ, ϕ)

θ̂ = train(
	  θ̂, θ_train, θ_val, Z_train, Z_val,
	  batchsize = batchsize, epochs = epochs, optimiser = ADAM(1e-5),
	  savepath = joinpath("intermediates", "Pretraining", "NotPretrained", "runs_N30"),
	  stopping_epochs = stopping_epochs
)
