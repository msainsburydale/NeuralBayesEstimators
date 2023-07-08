#!/bin/bash

set -e

echo "Do you wish to use a very low number of parameter configurations and epochs to quickly establish that the code is working (note that the generated results and plots will not exactly match those in the manuscript if you reply 'y')? (y/n) "
read quick_str

if [[ $quick_str == "y" ||  $quick_str == "Y" ]]; then
    quick=--quick
    ML=""
elif [[ $quick_str == "n" ||  $quick_str == "N" ]]; then
    quick=""
    ML=--ML
else
    echo "Please re-run and type y or n"
    exit 1
fi

for model in GaussianProcess/nuVaried GaussianProcess/nuFixed Schlather ConditionalExtremes
do
    echo ""
    echo "######## Starting simulation study for $model model ############"
    echo ""

    if [[ $model == "GaussianProcess/nuFixed" ]]; then
        ## Generate and plot the parameter configurations
        Rscript src/$model/Parameters.R $quick
        Rscript src/$model/PlotParameters.R
    fi

    ## Visualise field realisations
    julia --threads=auto --project=. src/VisualiseRealisations.jl --model=$model
    Rscript src/VisualiseRealisations.R --model=$model

    ## Train the neural estimators
    julia --threads=auto --project=. src/Train.jl --model=$model $quick --m="[1]"
    julia --threads=auto --project=. src/Train.jl --model=$model $quick --m="[1,10,30,75,150]" --deep

    ## Plot the loss functions
    Rscript src/PlotLoss.R --path=$model

    # Estimation
    if [[ $model == "ConditionalExtremes" ]]; then
        julia --threads=auto --project=. src/Estimate.jl --model=$model
    else
        julia --threads=auto --project=. src/Estimate.jl --model=$model $ML
    fi

    ## Results (plots, etc.)
    Rscript src/$model/Results.R
done
