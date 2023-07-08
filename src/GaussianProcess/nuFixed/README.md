# Gaussian process with fixed smoothness

In this study, the parameter configurations are generated in R (using code provided by Gerber and Nychka, 2021) and then loaded into Julia. Hence, `GaussianProcess/nuFixed` contains the scripts:

- `Parameters.R`: Generates parameter configurations and objects needed for data simulation (e.g., Cholesky factors). This file accepts a command-line argument `--quick`, a flag which, if present, drastically reduces the number of parameter configurations in each set.
- `Parameters.jl`: Defines functions for loading the above parameter configurations and data-simulation objects into julia, and defines a struct for storing these objects within julia.
