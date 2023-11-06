#!/bin/bash
unset R_HOME

set -e

echo "Do you wish to use a very low number of parameter configurations and epochs to quickly establish that the code is working? (y/n) "
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
julia --threads=auto --project=. src/Univariate/Train.jl --model=$model $quick --aggregation=logsumexp --m=10
julia --threads=auto --project=. src/Univariate/Estimate.jl --model=$model
Rscript src/Univariate/ResultsFixedSampleSize.R --model=$model

model=Normalsigma

echo ""
echo "######## Starting experiment on the effect of variable sample size: $model distribution ############"
echo ""

# Train the estimators
julia --threads=auto --project=. src/Univariate/Train.jl --model=$model $quick --m=5
julia --threads=auto --project=. src/Univariate/Train.jl --model=$model $quick --m=150
julia --threads=auto --project=. src/Univariate/Train.jl --model=$model $quick --m=-150

# Estimate and plot results
julia --threads=auto --project=. src/Univariate/Estimate.jl --model=$model
Rscript src/Univariate/ResultsVariableSampleSize.R --model=$model
