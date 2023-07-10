#!/bin/bash

set -e

echo ""
echo "######## Setting up ############"
echo ""

# Check that the data is installed
if [[ ! -f data/RedSea/redseatemperature.rdata ]]
then
    echo "The Red Sea data set has not been downloaded, or is in the wrong location. Please see the README for download instructions."
    exit 1
fi

# Install R dependencies
Rscript Dependencies.R


# Should we do a "quick" run through to check that things are working?
echo "Do you wish to use a very low number of parameter configurations and epochs to quickly establish that the code is working (note that the generated results and plots will not exactly match those in the manuscript if you reply 'y')? (y/n) "
read quick_str

if ! [[ $quick_str == "y" ||  $quick_str == "Y" || $quick_str == "n" ||  $quick_str == "N" ]]; then
    echo "Please re-run and type y or n"
    exit 1
fi

echo ""
echo "######## Starting the experiments... ############"
echo ""

# Each .sh files asks the user if quick = y/n. To automate this script,
# we pipe the above response to each .sh file
yes $quick_str | bash sh/Univariate.sh           # Section 2 and Section S2 of the Supplementary Material
yes $quick_str | bash sh/SimulationStudies.sh    # Section 3 and Section S7 of the Supplementary Material
yes $quick_str | bash sh/RedSea.sh               # Section 4
yes $quick_str | bash sh/SimulationOnTheFly.sh   # Section S4 of the Supplementary Material
yes $quick_str | bash sh/Pretraining.sh          # Section S5 of the Supplementary Material

echo ""
echo "######## Everything finished! ############"
echo ""
