# Source code

The methodology described in the manuscript has been developed into a user-friendly and well-documented `Julia` package, [NeuralEstimators.jl](https://github.com/msainsburydale/NeuralEstimators.jl), with an accompanying [`R` interface](https://github.com/msainsburydale/NeuralEstimators). The code in this repository is made available primarily for reproducibility purposes, and we encourage readers seeking to implement neural Bayes estimators to explore the package and its documentation.

The primary models considered in the manuscript are those discussed in Section 3, namely, the Gaussian process model, Schlather's max-stable model, and the spatial conditional extremes model; these models are represented by the folders `GaussianProcess`, `Schlather`, and `ConditionalExtremes`, respectively.  Each folder contains the scripts:
- `Parameters.jl`: Defines the prior measure, Ω(⋅), and intermediate objects needed for data simulation (e.g., distance matrices).
- `Simulation.jl`: Model-specific code for data simulation, which implicitly defines the assumed statistical model.
- `ML.jl`: Likelihood-based estimators (e.g., MAP and pairwise MAP estimators).
- `Results.R`: Generates figures and results.

The following scripts are common to the simulation studies of Section 3, and accept the command-line arguments listed below.
- `VisualiseRealisations.{jl/R}`: Simulates and plots realisations from the statistical model for a range of parameter configurations.
  - `--model`: A relative path to one of the simulation-study folders (e.g., `--model=Schlather`)
- `Train.jl`: Trains a neural estimator for a given statistical model.
  - `--model`
  - `--quick`: Flag controlling whether or not a computationally inexpensive run should be done.
  - `--m`: The sample size to use during training. If multiple samples sizes are given as a vector, multiple neural estimators will be trained.
  - `--deep`: Flag controlling whether to use a DeepSets neural estimator, or a neural estimator that is simply an average of one-at-a-time estimators.
- `Estimate.jl`: Parameter estimation over the test set, used to assess the estimators.
  - `--model`
  - `--quick`
  - `--ML`: Flag controlling whether or not likelihood estimation should be performed.

Other studies in the manuscript follow a similar pattern. Note that the code for the Gaussian process model with fixed smoothness handles the parameter configurations slightly differently; see [`GaussianProcess/nuFixed/README`](https://github.com/msainsburydale/NeuralBayesEstimators/tree/master/src/GaussianProcess/nuFixed) for details.

<!-- Two additional scripts that are common to all experiments are:

- `PlotLoss.R`: Plots the loss-function evolution during training. Takes the command-line argument `--path`, which should point to a folder containing an arbitrary number of subfolders named `runs_*`; the loss during training for each of these folders will be plotted in a single figure.
- `Architecture.jl`: The common architecture used for the neural estimators. -->


# Controlling shell scripts

The experiments are conducted in separate Julia sessions and coordinated using `.sh` files (in the folder `sh`). This prevents garbage-collection issues and dependencies that could affect the comparison studies.
