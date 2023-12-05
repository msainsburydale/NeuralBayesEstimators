include(joinpath(pwd(), "src", "ConditionalExtremes", "Simulation.jl"))


function reshapeZ(Z::A, ξ) where {T <: Number, A <: AbstractArray{T, 2}}
    # Just need to add a singleton second dimension
    return reshape(Z, size(Z, 1), 1, size(Z, 2))
end
