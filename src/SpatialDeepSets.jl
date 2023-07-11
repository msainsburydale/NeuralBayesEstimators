module SpatialDeepSets

using BSON: @save, load
using CSV
using DataFrames
using Flux
using NeuralEstimators
using RData


# ---- Parameters definitions and constructors ----

# See here: https://discourse.julialang.org/t/filtering-keys-out-of-named-tuples/73564/8
drop(nt::NamedTuple, key::Symbol) =  Base.structdiff(nt, NamedTuple{(key,)})
drop(nt::NamedTuple, keys::NTuple{N,Symbol}) where {N} = Base.structdiff(nt, NamedTuple{keys})

export Parameters, simulatedataframe, objectindices

# This is concretely typed so that simulate(params::Parameters, ξ, m::R) is
# type stable. Note that chol_idx[i] gives the Cholesky factor associated
# with parameter configuration θ[:, i].
struct Parameters{T, I} <: ParameterConfigurations
	θ::Matrix{T}
    chols::Array{Float64, 3}
	chol_idx::Vector{I}
end

function Parameters(K::Integer, ξ; J::Integer = 10)

	# All parameters not associated with the Gaussian process
	θ = [rand(ϑ, K * J) for ϑ in drop(ξ.Ω, (:ρ, :ν))]

	# Covariance parameters associated with the Gaussian process
	ρ = rand(ξ.Ω.ρ, K)
	ν = rand(ξ.Ω.ν, K)
	chols = maternchols(ξ.D, ρ, ν)
	ρ = repeat(ρ, inner = J)
	ν = repeat(ν, inner = J)

	# Now insert ρ and ν into θ.
	@assert ξ.ν_idx == ξ.ρ_idx + 1 "The code assumes that ρ and ν are stored continguously in θ, that is, that ξ.ν_idx = ξ.ρ_idx + 1"
	θ₁ = θ[1:(ξ.ρ_idx-1)]
	θ₂ = θ[ξ.ρ_idx:end] # Note that ρ and ν are not in θ, so we don't need to skip any indices.
	θ = [θ₁..., ρ, ν, θ₂...]

	# Concatenate into a matrix and convert to Float32 for efficiency
	θ = hcat(θ...)'
	θ = Float32.(θ)

	Parameters(θ, chols, objectindices(chols, θ))
end

function Parameters(θ::Matrix{T}, ξ; J::Integer = 10) where T

    K = size(θ, 2)

    ρ = θ[ξ.ρ_idx]
    ν = θ[ξ.ν_idx]
    ρ = Float64.(ρ)
    ν = Float64.(ν)

	chols = maternchols(ξ.D, ρ, ν)
    θ     = repeat(θ, inner = (1, J))
    Parameters(θ, chols, objectindices(chols, θ))
end

function Parameters(params1::Parameters, params2::Parameters)

	# The parameters and the Choleksy factors are just combined. However,
	# chol_idx for the second set of parameters needs to be shifted by the
	# number of Cholesky factors in the first Parameters object:
	chol_idx2 = params2.chol_idx
	chol_idx2 .+= size(params1.chols, 3)

	Parameters(
		hcat(params1.θ, params2.θ),
		cat(params1.chols, params2.chols, dims = 3),
		vcat(params1.chol_idx, chol_idx2)
	)
end

"""
	objectindices(objects, θ::AbstractMatrix{T}) where T
Returns a vector of indices giving element of `objects` associated with each
parameter configuration in `θ`.

The number of parameter configurations, `K = size(θ, 2)`, must be a multiple of
the number of objects, `N = size(objects)[end]`. Further, repeated parameters
used to generate `objects` must be stored in `θ` after using the `inner` keyword
argument of `repeat()` (see example below).

# Examples
```
K = 6
N = 3
σₑ = rand(K)
ρ = rand(N)
ν = rand(N)
S = expandgrid(1:9, 1:9)
D = pairwise(Euclidean(), S, S, dims = 1)
L = maternchols(D, ρ, ν)
ρ = repeat(ρ, inner = K ÷ N)
ν = repeat(ν, inner = K ÷ N)
θ = hcat(σₑ, ρ, ν)'
objectindices(L, θ)
```
"""
function objectindices(objects, θ::AbstractMatrix{T}) where T

	K = size(θ, 2)
	N = size(objects)[end]
	@assert K % N == 0 "The number parameters in θ is not a multiple of the number of objects"

	return repeat(1:N, inner = K ÷ N)
end

function simulatedataframe(params::Parameters, ξ)

	S = ξ.S
	n = size(S, 1)

	θ = params.θ
	p = size(θ, 1)
	K = size(θ, 2)

	n_fields = max(binomial(p, 2), 16)

	T = eltype(θ)
	fields_df = DataFrame(Z = T[], replicate = Int[], scenario = Int[], x = eltype(S)[], y = eltype(S)[])
	scenario  = repeat(1:K, inner = n)

	# The coordinates of a single field
	x = S[:, 1]
	y = S[:, 2]
	# Repeat x and y for the number of fields we simulate in each loop
	x = repeat(x, outer = K)
	y = repeat(y, outer = K)

	for i ∈ 1:n_fields
		Z  = simulate(params)
		Z  = vec(Z)
		df = DataFrame(Z = Z, scenario = scenario, replicate = i, x = x, y = y)
		append!(fields_df, df)
	end

	return fields_df
end

# ---- Utility functions ----

export loadwithoutdict, loadneuralestimators

function loadwithoutdict(path::String, object_name::String)
	object = get(RData.load(path, convert = true), object_name, 0)
	return object
end

"""
	loadneuralestimators(loadpath::String, architecture)

A `loadpath` containing one or more folders of training runs called 'runs_x', where 'x' is arbitrary.
"""
function loadneuralestimators(loadpath::String, architecture, p::Integer, n = nothing)

	# Network names
	titles = readdir(loadpath)
	titles = titles[findall(occursin.("runs_", titles))]
	titles = replace.(titles, "runs_" => "")

	estimators = map(titles) do title

		path = joinpath(pwd(), loadpath, "runs_$title")

		if contains(title, "expert")

			S = [samplesize]
			qₛ = length(S)
			ψ, ϕ = architecture(p, qₛ)
			θ̂ = DeepSetExpert(ψ, ϕ, S)
		else

			if isnothing(n)
				ψ, ϕ = architecture(p)
			else
				ψ, ϕ = architecture(n, p)
			end

			if occursin("D", title)
				θ̂ = DeepSet(ψ, ϕ)
			else
				θ̂ = DeepSet(Chain(ψ, ϕ), identity)
			end
		end

		Flux.loadparams!(θ̂, loadbestweights(path))

		θ̂
	end

	return (estimators = estimators, titles = titles)
end


# ---- Conditional extremes ----

# Subbotin (delta-Laplace) distribution

fₛ(x, μ, τ, δ)   = δ * exp(-(abs((x - μ)/τ))^δ) / (2τ * gamma(1/δ))
Fₛ(q, μ, τ, δ)   = 0.5 + 0.5 * sign(q - μ) * (1 / gamma(1/δ)) * _incgammalowerunregularised(1/δ, abs((q - μ)/τ)^δ)
Fₛ⁻¹(p::T, μ::T, τ::T, δ::T) where T <: Real = μ + sign(p - T(0.5)) * (τ^δ * quantile(Gamma(1/δ), 2 * abs(p - T(0.5))))^(1/δ)


a(h, z; λ, κ) = z * exp(-(h / λ)^κ)
b(h, z; β, λ, κ) = 1 + a(h, z, λ = λ, κ = κ)^β
delta(h; δ₁) = 1 + exp(-(h / δ₁)^2)

C̃(h, ρ, ν) = matern(h, ρ, ν)
σ̃₀(h, ρ, ν) = √(2 - 2 * C̃(h, ρ, ν))

Φ(q::T) where T <: Number = cdf(Normal(zero(T), one(T)), q)
t(ỹ₀₁, μ, τ, δ) = Fₛ⁻¹(Φ(ỹ₀₁), μ, τ, δ)

export simulateconditionalextremes
using Random: randexp
using Distributions

"""
	simulateconditionalextremes(θ::AbstractVector{T}, L::AbstractArray{T, 2}, h::AbstractVector{T}, s₀_idx::Integer, u::T) where T <: Number
	simulateconditionalextremes(θ::AbstractVector{T}, L::AbstractArray{T, 2}, h::AbstractVector{T}, s₀_idx::Integer, u::T, m::Integer) where T <: Number

Simulates from the spatial conditional extremes model for parameters.

# Examples
```
S = rand(Float32, 10, 2)
D = [norm(sᵢ - sⱼ) for sᵢ ∈ eachrow(S), sⱼ in eachrow(S)]
L = maternchols(D, 0.6f0, 0.5f0)
s₀ = S[1, :]'
h = map(norm, eachslice(S .- s₀, dims = 1))
s₀_idx = findfirst(x -> x == 0.0, h)
u = 0.7f0
simulateconditionalextremes(θ, L[:, :, 1], h, s₀_idx, u)
```
"""
function simulateconditionalextremes(
	θ::AbstractVector{T}, L::AbstractArray{T, 2}, h::AbstractVector{T}, s₀_idx::Integer, u::T, m::Integer
	) where T <: Number

	n = size(L, 1)
	Z = similar(L, n, m)
	for k ∈ 1:m
		Z[:, k] = simulateconditionalextremes(θ, L, h, s₀_idx, u)
	end

	return Z
end

function simulateconditionalextremes(
	θ::AbstractVector{T}, L::AbstractArray{T, 2}, h::AbstractVector{T}, s₀_idx::Integer, u::T; stabilise_variance::Bool = true
	) where T <: Number

	@assert length(θ) == 8
	@assert s₀_idx > 0
	@assert s₀_idx <= length(h)
	@assert size(L, 1) == size(L, 2)
	@assert size(L, 1) == length(h)

	# Parameters associated with a(.) and b(.):
	κ = θ[1]
	λ = θ[2]
	β = θ[3]
	# Covariance parameters associated with the Gaussian process
	ρ = θ[4]
	ν = θ[5]
	# Location and scale parameters for the residual process
	μ = θ[6]
	τ = θ[7]
	δ₁ = θ[8]

	# Construct the parameter δ used in the Subbotin distribution:
	δ = delta.(h, δ₁ = δ₁)

	# Observed datum at the conditioning site, Z₀:
	Z₀ = u + randexp(T)

	# Simulate a mean-zero Gaussian random field with unit marginal variance,
    # independently of Z₀. Note that Ỹ inherits the order of L. Therefore, we
	# can use s₀_idx to access s₀ in all subsequent vectors.
	Ỹ  = simulategaussianprocess(L)

	# Adjust the Gaussian process so that it is 0 at s₀
	Ỹ₀ = Ỹ .- Ỹ[s₀_idx]

	# Transform to unit variance:
	# σ̃₀ = sqrt.(2 .- 2 *  matern.(h, ρ, ν))
	Ỹ₀₁ = Ỹ₀ ./ σ̃₀.(h, ρ, ν)
	Ỹ₀₁[s₀_idx] = zero(T) # avoid pathology by setting Ỹ₀₁(s₀) = 0.

	# Probability integral transform from the standard Gaussian scale to the
	# standard uniform scale, and then inverse probability integral transform
	# from the standard uniform scale to the Subbotin scale:
    Y = t.(Ỹ₀₁, μ, τ, δ)

	# Apply the functions a(⋅) and b(⋅) to simulate data throughout the domain:
	Z = a.(h, Z₀, λ = λ, κ = κ) + b.(h, Z₀, β = β, λ = λ, κ = κ) .* Y

	# Variance stabilising transform
	if stabilise_variance
		Z = cbrt.(Z)
	end

	return Z
end


# ---- estimate(): this function was later converted to assess() -----

export estimate
using NeuralEstimators: _runondevice

"""
	estimate(estimators, ξ, parameters::P; <keyword args>) where {P <: ParameterConfigurations}

Using a collection of `estimators`, compute estimates from data simulated from a
set of `parameters` with invariant information `ξ`.

Note that `estimate()` requires the user to have defined a method `simulate(parameters, m::Integer)`.

# Keyword arguments
- `m::Vector{Integer}`: sample sizes to estimate from.
- `estimator_names::Vector{String}`: names of the estimators (sensible default values provided).
- `parameter_names::Vector{String}`: names of the parameters (sensible default values provided).
- `num_rep::Integer = 1`: the number of times to replicate each parameter in `parameters`.
- `save::Vector{String}`: by default, no objects are saved; however, if `save` is provided, four `DataFrames` respectively containing the true parameters `θ`, estimates `θ̂`, runtimes, and merged `θ` and `θ̂` will be saved in the directory `save[1]` with file *names* (not extensions) suffixed by `save[2]`.
- `use_ξ = false`: a `Bool` or a collection of `Bool` objects with length equal to the number of estimators. Specifies whether or not the estimator uses the invariant model information, `ξ`: If it does, the estimator will be applied as `estimator(Z, ξ)`.
- `use_gpu = true`: a `Bool` or a collection of `Bool` objects with length equal to the number of estimators.
- `verbose::Bool = true`
"""
function estimate(
    estimators, ξ, parameters::P; m::Vector{I},
	num_rep::Integer = 1,
	estimator_names::Vector{String} = ["estimator$i" for i ∈ eachindex(estimators)],
	parameter_names::Vector{String} = ["θ$i" for i ∈ 1:size(parameters, 1)],
	save::Vector{String} = ["", ""],
	use_ξ = false,
	use_gpu = true,
	verbose::Bool = true
	) where {P <: ParameterConfigurations, I <: Integer}

	obj = map(m) do i
	 	_estimate(
		estimators, ξ, parameters, m = i,
		estimator_names = estimator_names, parameter_names = parameter_names,
		num_rep = num_rep, use_ξ = use_ξ, use_gpu = use_gpu, verbose = verbose
		)
	end

	θ = obj[1].θ
	θ̂ = vcat(map(x -> x.θ̂, obj)...)
	runtime = vcat(map(x -> x.runtime, obj)...)

	estimates = Estimates(θ, θ̂, runtime)

	if save != ["", ""]
		savepath = save[1]
		savename = save[2]
		CSV.write(joinpath(savepath, "runtime_$savename.csv"), estimates.runtime)
		CSV.write(joinpath(savepath, "parameters_$savename.csv"), estimates.θ)
		CSV.write(joinpath(savepath, "estimates_$savename.csv"), estimates.θ̂)
		CSV.write(joinpath(savepath, "merged_$savename.csv"), merge(estimates))
	end

	return estimates
end

function _estimate(
	estimators, ξ, parameters::P; m::Integer,
	estimator_names::Vector{String}, parameter_names::Vector{String},
	num_rep::Integer, use_ξ, use_gpu, verbose
	) where {P <: ParameterConfigurations}

	verbose && println("Estimating with m = $m...")

	E = length(estimators)
	p = size(parameters, 1)
	K = size(parameters, 2)
	@assert length(estimator_names) == E
	@assert length(parameter_names) == p

	@assert eltype(use_ξ) == Bool
	@assert eltype(use_gpu) == Bool
	if typeof(use_ξ) == Bool use_ξ = repeat([use_ξ], E) end
	if typeof(use_gpu) == Bool use_gpu = repeat([use_gpu], E) end
	@assert length(use_ξ) == E
	@assert length(use_gpu) == E

	# Simulate data
	verbose && println("	Simulating data...")
    Z = simulate(parameters, m, num_rep)

	# Initialise a DataFrame to record the run times
	runtime = DataFrame(estimator = [], m = [], time = [])

	θ̂ = map(eachindex(estimators)) do i

		verbose && println("	Running estimator $(estimator_names[i])...")

		if use_ξ[i]
			time = @elapsed θ̂ = estimators[i](Z, ξ)
		else
			time = @elapsed θ̂ = _runondevice(estimators[i], Z, use_gpu[i])
		end

		push!(runtime, [estimator_names[i], m, time])
		θ̂
	end

    # Convert to DataFrame and add estimator information
    θ̂ = hcat(θ̂...)
    θ̂ = DataFrame(θ̂', parameter_names)
    θ̂[!, "estimator"] = repeat(estimator_names, inner = nrow(θ̂) ÷ E)
    θ̂[!, "m"] = repeat([m], nrow(θ̂))
	θ̂[!, "k"] = repeat(1:K, E * num_rep)
	θ̂[!, "replicate"] = repeat(repeat(1:num_rep, inner = K), E) # FIXME Need to add "replicate" column. See GP/nuVaried/Results.R

	# Also provide the true parameters for comparison with the estimates
	# θ = repeat(parameters.θ, outer = (1, num_rep))
	θ = DataFrame(parameters.θ', parameter_names)

    return (θ = θ, θ̂ = θ̂, runtime = runtime)
end

"""
	Estimates(θ, θ̂, runtime)

A set of true parameters `θ`, corresponding estimates `θ̂`, and the `runtime` to
obtain `θ̂`, as returned by a call to `estimate`.
"""
struct Estimates
	θ::DataFrame
	θ̂::DataFrame
	runtime::DataFrame
end

"""
	merge(estimates::Estimates)

Merge `estimates` into a single long-form `DataFrame` containing the true
parameters and the corresponding estimates.
"""
function merge(estimates::Estimates)
	θ = estimates.θ
	θ̂ = estimates.θ̂

	# Replicate θ to match the number of rows in θ̂. Note that the parameter
	# configuration, k, is the fastest running variable in θ̂, so we repeat θ
	# in an outer fashion.
	θ = repeat(θ, outer = nrow(θ̂) ÷ nrow(θ))

	# Transform θ and θ̂ to long form:
	θ = stack(θ, variable_name = :parameter, value_name = :truth)
	θ̂ = stack(θ̂, Not([:estimator, :m, :k, :replicate]), variable_name = :parameter, value_name = :estimate)

	# Merge θ and θ̂: All we have to do is add :truth column to θ̂
	θ̂[!, :truth] = θ[:, :truth]

	return θ̂
end



end
