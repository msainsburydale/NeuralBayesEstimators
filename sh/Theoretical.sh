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

model=Uniform

echo ""
echo "######## Starting experiment on approximating the Bayes estimator: $model distribution ############"
echo ""

# Show that NN estimators can approximate the Bayes estimator for a one-parameter model
julia --threads=auto --project=. src/Theoretical/Train.jl --model=$model $quick --aggregation=logsumexp --qt=100 --m=10
julia --threads=auto --project=. src/Theoretical/Estimate.jl --model=$model
Rscript src/Theoretical/Results.R --model=$model
Rscript src/Theoretical/ResultsFixedSetSize.R --model=$model
