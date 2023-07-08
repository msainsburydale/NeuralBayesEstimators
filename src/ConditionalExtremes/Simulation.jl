using NeuralEstimators
import NeuralEstimators: simulate
using SpatialDeepSets
using SpatialDeepSets: fₛ, Fₛ, Fₛ⁻¹
using Folds

"""
	simulate(parameters::Parameters)
	simulate(parameters::Parameters, m::Integer)
	simulate(parameters::Parameters, m::R) where {R <: AbstractRange{I}} where I <: Integer

Simulates `m` fields from a spatial conditional extremes model for each of the
given covariance `parameters`.

If `m` is not provided, a single field is simulated for each parameter
configuration, and the return type is an array with the last dimension
corresponding to the parameters. If `m` is provided, `m` fields
are simulated for each parameter configuration, and the return type is a vector
of arrays equal in length to the number of parameter configurations, and with
the fourth dimension of the array containing the field replicates.
"""
function simulate(params::Parameters, m::R) where {R <: AbstractRange{I}} where I <: Integer

	K = size(params, 2)
	m̃ = rand(m, K)

	θ         = Float64.(params.θ)
	chols     = params.chols
	chol_idx  = params.chol_idx

	# NB Folds.map() is type unstable. I've open an issue with the package, but
	# have not received a response yet. To improve efficiency, I may need to use
	# an alternative parallel mapping function.
	Z = Folds.map(1:K) do k
		L = view(chols, :, :, chol_idx[k])
		z = simulateconditionalextremes(θ[:, k], L, ξ.h, ξ.s₀_idx, ξ.u, m̃[k]) # nxm matrix
		z = reshapeZ(z, ξ)
		z = Float32.(z)
		z
	end

    return Z
end

simulate(parameters::Parameters, m::Integer) = simulate(parameters, range(m, m))
simulate(parameters::Parameters) = stackarrays(simulate(parameters, 1))

"""
	reshapeZ(Z::A, ξ) where {T <: Number, A <: AbstractArray{T, 2}}
Reshape the data to the appropriate form for the given architecture and
application. Note that this method may be overwritten in other scripts.
"""
function reshapeZ(Z::A, ξ) where {T <: Number, A <: AbstractArray{T, 2}}
	n = size(Z, 1)
	m = size(Z, 2)
	@assert isqrt(n) == √n "n should be a whole number"
	return reshape(Z, isqrt(n), isqrt(n), 1, m)
end
