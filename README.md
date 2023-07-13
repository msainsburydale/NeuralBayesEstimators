# Source code for "Likelihood-Free Parameter Estimation with Neural Bayes Estimators"

This repository contains code for reproducing the results in ["Likelihood-Free Parameter Estimation with Neural Bayes Estimators" (Sainsbury-Dale, Zammit-Mangion, and Huser, 2022)](https://arxiv.org/abs/2208.12942).

The methodology described in the manuscript has been developed into a user-friendly and well-documented Julia package, [NeuralEstimators.jl](https://github.com/msainsburydale/NeuralEstimators.jl), with an accompanying [R interface](https://github.com/msainsburydale/NeuralEstimators). The code in this repository is made available primarily for reproducibility purposes, and we encourage readers seeking to implement neural Bayes estimators to explore the package and its documentation.  

## Repository structure

We first briefly describe the repository structure, although an understanding of this structure is not needed for reproducing the results. The repository is organised into folders containing source code (`src`), intermediate objects generated from the source code (`intermediates`), figures (`img`), results (`results`), and controlling shell scripts that weave everything together (`sh`). Each folder is further divided into the following tree structure, where each branch is associated with one component of the manuscript:

```bash
├── Univariate          (Section 2)
├── GaussianProcess
│   ├── nuFixed         (Section S7 of the Supplementary Material)
│   ├── nuVaried        (Section 3.2)
├── Schlather           (Section 3.3)
├── ConditionalExtremes (Section 3.4)
├── RedSea              (Section 4)
├── SimulationOnTheFly  (Section S4 of the Supplementary Material)
```

Further details are given in [`src/README.md`](https://github.com/msainsburydale/NeuralBayesEstimators/tree/master/src).

## Instructions

First, download this repository and navigate to its top-level directory within the command line (i.e., `cd` to wherever you installed the repository).

### Data

The Red Sea data set, analysed in Section 4, was too large (a few hundred Mb in total) to be stored on GitHub. To automatically download it and place it into the correct location, run `bash sh/data.sh`. If the data are not downloading as expected, please download them from [here](https://zenodo.org/record/8134200) and place them into the folder `data/RedSea`.
<!-- [here](https://hpc.niasra.uow.edu.au/ckan/dataset/red_sea_temperature) (to download the file, click "Explore" > "Go to resource") and place them into the folder `data/RedSea`. -->


Note that checks at the beginning of the replication script will immediately notify the user if these data are not present.

### Software dependencies

Before installing the software dependencies, users may wish to setup a `conda` environment, so that the dependencies of this repository do not affect the user's current installation. To create a `conda` environment, run the following command in terminal:

```
conda create -n NeuralBayesEstimators -c conda-forge julia=1.7.1 r-base nlopt
```

Then activate the `conda` environment with:

```
conda activate NeuralBayesEstimators
```

The above `conda` environment installs Julia and R automatically. If you do not wish to use a `conda` environment, you will need to install Julia and R manually if they are not already on your system:  

- Install Julia 1.7.1. (See [here](https://julialang.org/downloads/).)
  - Ensure that your system can find the `julia` executable (this usually needs to be done manually; see, e.g., [here](https://julialang.org/downloads/platform/#linux_and_freebsd)) by entering `julia` in terminal, which should open the Julia REPL (run `exit()` to leave the REPL).
- Install R >= 4.0.0. (See [here](https://www.r-project.org/).)

Once Julia and R are setup, install package dependencies as follows:

- In terminal, navigate (i.e., `cd`) to the top level of this repository, and enter:
  - `julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'`. This will install all Julia package dependencies given in `Project.toml` and `Manifest.toml`.
  - `Rscript Dependencies.R`. This will install all R package dependencies given in `Dependencies.txt`.


### Hardware requirements

The fast construction of neural Bayes estimators requires graphical processing units (GPUs). Hence, although the code in this repository will run without a GPU (i.e., it will run on the CPU), we recommend that the user run this code on a workstation with a GPU. Note that running the "quick" version of the code (see below) is still fast even on the CPU, so the code can still be tested without a GPU.

The Red-Sea application study of Section 4 is memory intensive; at least 128GB of of CPU RAM (or RAM + swap) is needed to run the full (i.e., non "quick") version of this study. If this is an issue, please comment out the line containing `bash sh/RedSea.sh` in the replication script, `sh/all.sh`. (Comments in `.sh` files are made with `#`.)

### Reproducing the results

The replication script is `sh/all.sh`, invoked using `bash sh/all.sh` from the top level of this repository. For all studies, the replication script will automatically train the neural estimators, generate estimates from both the neural and likelihood-based estimators, and populate the `img` and `results` folders with the figures and results of the manuscript.

The nature of our experiments means that the run time for reproducing the results of the manuscript is substantial (2-3 days in total). When running the replication script, the user will be prompted with an option to quickly establish that the code is working by using a small number of parameter configurations and epochs. Our envisioned workflow is to establish that the code is working with this "quick" option, clear the populated folders by simply entering `bash sh/clear.sh`, and then run the code in full (possibly over the weekend).

Note that the replication script is clearly presented and commented; hence, one may easily "comment out" sections to produce a subset of the results. (Comments in `.sh` files are made with `#`.)

#### Minor reproducibility difficulties

When training neural networks on the GPU, there is some some unavoidable non-determinism: see [here](https://discourse.julialang.org/t/flux-reproducibility-of-gpu-experiments/62092). This does not significantly affect the "story" of the final results, but there may be some slight differences each time the code is executed.
