# To reduce code repetition, I have used a single DataSimulation.jl file
# for both ν fixed and ν unknown. Here, we simply source that file.
# Note that we still include a DataSimulation.jl file within each
# subfolder for consistency with the assumed structure of the repo.
include(joinpath(pwd(), "src", "GaussianProcess", "Simulation.jl"))