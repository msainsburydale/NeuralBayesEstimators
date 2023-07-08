#!/bin/bash

set -e

echo "Do you wish to use a very low number of parameter configurations and epochs to quickly establish that the code is working (note that the generated results and plots will not exactly match those in the manuscript if you reply 'y')? (y/n) "
read quick_str
if [[ $quick_str == "y" ||  $quick_str == "Y" ]]; then
    quick=--quick
elif [[ $quick_str == "n" ||  $quick_str == "N" ]]; then
    quick=""
else
    echo "Please re-run and type y or n"
    exit 1
fi

model=Normalsigma
echo ""
echo "######## Starting experiment on the effect of variable sample size: $model distribution ############"
echo ""

# Train the estimators
qt=64
julia --threads=auto --project=. src/Theoretical/Train.jl --model=$model $quick --qt=$qt --prior=informative --m=-150 --aggregation=mean
julia --threads=auto --project=. src/Theoretical/Train.jl --model=$model $quick --qt=$qt --prior=informative --m=5    --aggregation=mean
julia --threads=auto --project=. src/Theoretical/Train.jl --model=$model $quick --qt=$qt --prior=informative --m=150  --aggregation=mean

# Estimate and plot results
julia --threads=auto --project=. src/Theoretical/Estimate.jl --model=$model --prior=$prior
Rscript src/Theoretical/ResultsVariableSetSize.R --model=$model
