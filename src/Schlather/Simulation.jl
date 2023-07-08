using NeuralEstimators
import NeuralEstimators: simulate
using SpatialDeepSets
using Folds

function simulate(parameters::Parameters, m::R) where {R <: AbstractRange{I}} where I <: Integer

	P = size(parameters, 2)
	m̃ = rand(m, P)

	chols    = parameters.chols
	chol_idx = parameters.chol_idx

	# NB Folds.map() is type unstable. I've open an issue with the package, but
	# have not received a response yet. To improve efficiency, I may need to use
	# an alternative parallel mapping function.
	Z = Folds.map(1:P) do i
		L = view(chols, :, :, chol_idx[i])
		z = simulateschlather(L, m̃[i])
		z = Float32.(z)
		z
	end
	n = size(chols, 1)
	Z = reshape.(Z, isqrt(n), isqrt(n), 1, :) # assumes a square domain
	return Z
end

simulate(parameters::Parameters, m::Integer) = simulate(parameters, range(m, m))
simulate(parameters::Parameters) = stackarrays(simulate(parameters, 1))
