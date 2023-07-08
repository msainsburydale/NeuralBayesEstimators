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

echo ""
echo "######## Starting experiment on the effect of 'on-the-fly' simulation ############"
echo ""

julia --threads=auto --project=. src/SimulationOnTheFly/Train.jl $quick
Rscript src/SimulationOnTheFly/Results.R
