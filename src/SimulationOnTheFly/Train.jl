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
using SpatialDeepSets
using DataFrames
using CSV
using Random: seed!
using Flux

model = "GaussianProcess/nuFixed"
include(joinpath(pwd(), "src/$model/Parameters.jl"))
include(joinpath(pwd(), "src/$model/Simulation.jl"))

epochs = quick ? 10 : 150

# ---- Define the prior distribution ----

# Use the same prior as we do for the Gaussian process with unknown smoothness

Ω = (
	σ = Uniform(0.1, 1),
	ρ = Uniform(2, 10)
)
parameter_names = String.(collect(keys(Ω)))
S = expandgrid(1:16, 1:16)
D = pairwise(Euclidean(), S, S, dims = 1)
ξ = (
	Ω = Ω, S = S, D = D, p = length(Ω),
	parameter_names = parameter_names,
	ρ_idx = findfirst(parameter_names .== "ρ"),
	invtransform = identity
)

# ---- Construct the training and validation parameters ----

K_train = 10_000
K_val   = K_train ÷ 5
if quick
	K_train = K_train ÷ 100
	K_val   = K_val ÷ 100
end

@info "Sampling validation parameters..."
θ_val   = Parameters(K_val, ξ)
@info "Sampling training parameters..."
θ_train = Parameters(K_train, ξ)

# ---- Architectures to use in the study ----

# Network with an order-of-magnitude fewer parameters (~69,000) than the architecture used by Gerber and Nychka (2021)
function small_network(p)
	Chain(
		Conv((10, 10), 1 => 16, relu),
		Conv((5, 5), 16 => 32, relu),
		Conv((3, 3), 32 => 64, relu),
		Flux.flatten,
		Dense(64, 500, relu),
		Dense(500, p)
	)
end

# Architecture used by Gerber and Nychka (2021): Assumes a 16x16 matrix as input. 635_742 parameters when p = 2.
function large_network(p)
	Chain(
		Conv((10, 10), 1 => 128, relu),
		Conv((5, 5), 128 => 128, relu),
		Conv((3, 3), 128 => 128, relu),
		Flux.flatten,
		Dense(128, 500, relu),
		Dense(500, p)
	)
end

architectures = (smallNetwork = small_network, largeNetwork = large_network)

# ---- Training ----
m = 1
epochs_per_Z_refreshes = (everyEpoch = 1, someEpochs = 30, noEpochs = 10000)

for arch_name ∈ keys(architectures)
	for epochs_per_Z_refresh_name ∈ keys(epochs_per_Z_refreshes)

		epochs_per_Z_refresh = epochs_per_Z_refreshes[epochs_per_Z_refresh_name]

		@info "Training the $arch_name network, refreshing the training set every $epochs_per_Z_refresh epochs"

		seed!(1)
		θ̂ = architectures[arch_name](ξ.p)
		θ̂ = DeepSet(θ̂, identity)

		θ̂ = train(
		      θ̂, θ_train, θ_val, simulate, m = m,
		      epochs_per_Z_refresh = epochs_per_Z_refreshes[epochs_per_Z_refresh_name],
		      epochs = epochs, batchsize = 256, stopping_epochs = epochs, # no early stopping
		      savepath = "intermediates/SimulationOnTheFly/$arch_name/$epochs_per_Z_refresh_name/runs_N1"
		)
	end
end


# Fixed training data but larger amount of training data than validation
Z_train = simulate(θ_train, ξ, 3m)
Z_val   = simulate(θ_val, ξ, m)
for arch_name ∈ keys(architectures)
	@info "Training the $arch_name network with recycled training data"

	seed!(1)
	θ̂ = architectures[arch_name](ξ.p)
	θ̂ = DeepSet(θ̂, identity)

	θ̂ = train(
	      θ̂, θ_train, θ_val, Z_train, Z_val,
	      epochs = epochs, batchsize = 256, stopping_epochs = epochs, # no early stopping
	      savepath = "intermediates/SimulationOnTheFly/$arch_name/recycled/runs_N1"
	)
end
