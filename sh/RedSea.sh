#!/bin/bash
unset R_HOME

set -e

if [[ ! -f data/RedSea/redseatemperature.rdata ]]
then
    echo "The Red Sea data set has not been downloaded, or is in the wrong location. Please see the README for download instructions."
    exit 1
fi


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

echo ""
echo "######## Starting application study for the Red Sea data set ############"
echo ""

# Data processing and defining the set of parameter configurations
Rscript src/RedSea/DataProcessing.R --plot
Rscript src/RedSea/Parameters.R $quick

# Empirical Bootstrap
julia --threads=auto --project=. src/RedSea/EmpiricalEstimatorBootstrap.jl

# Train and estimate with each neural estimator described in the manuscript
julia --threads=auto --project=. src/RedSea/Train.jl $quick --arch=CNN --data_type=regular
Rscript src/PlotLoss.R --path=RedSea/CNN
julia --threads=auto --project=. src/RedSea/Estimate.jl --arch=CNN --data_type=regular
Rscript src/RedSea/Results.R --arch=CNN --data_type=regular

for datatype in regular irregular
do
    julia --threads=auto --project=. src/RedSea/Train.jl $quick --arch=DNN  --data_type=$datatype
    Rscript src/PlotLoss.R --path=RedSea/DNN$datatype
    julia --threads=auto --project=. src/RedSea/Estimate.jl --arch=DNN  --data_type=$datatype
    Rscript src/RedSea/Results.R --arch=DNN --data_type=$datatype
done

# Plots that involve all of the neural estimators
Rscript src/RedSea/ResultsCombined.R
