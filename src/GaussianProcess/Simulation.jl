using NeuralEstimators
import NeuralEstimators: simulate
using SpatialDeepSets
using Folds

"""
	simulate(parameters::Parameters)
	simulate(parameters::Parameters, m::Integer)
	simulate(parameters::Parameters, m::R) where {R <: AbstractRange{I}} where I <: Integer

Simulates `m` fields from a Gaussian process for each of the given covariance `parameters`.

If `m` is not provided, a single field is simulated for each parameter
configuration, and the return type is an array with the last dimension
corresponding to the parameters. If `m` is provided, `m` fields
are simulated for each parameter configuration, and the return type is a vector
of arrays equal in length to the number of parameter configurations, and with
the fourth dimension of the array containing the field replicates.

This function assumes that the nugget standard deviation parameters is stored in the first row of `parameters.θ`.
"""
function simulate(parameters::Parameters, m::R) where {R <: AbstractRange{I}} where I <: Integer

	P = size(parameters, 2)
	m̃ = rand(m, P)

	σ        = Float64.(parameters.θ[1, :])
	chols    = parameters.chols
	chol_idx = parameters.chol_idx

	# NB Folds.map() is type unstable. I've open an issue with the package, but
	# have not received a response yet. To improve efficiency, I may need to use
	# an alternative parallel mapping function.
	Z = Folds.map(1:P) do i
		L = view(chols, :, :, chol_idx[i])
		L = Matrix(L)
		z = simulategaussianprocess(L, m̃[i])
		z = z + σ[i] * randn(size(z)...) # add measurement error
		z = Float32.(z)
		z
	end
	n = size(chols, 1)
	Z = reshape.(Z, isqrt(n), isqrt(n), 1, :) # assumes a square domain
	return Z
end

simulate(parameters::Parameters, m::Integer) = simulate(parameters, range(m, m))
simulate(parameters::Parameters) = stackarrays(simulate(parameters, 1))
