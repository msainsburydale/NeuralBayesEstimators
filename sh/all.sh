#!/bin/bash
unset R_HOME

set -e

echo ""
echo "######## Setting up ############"
echo ""

# Should we do a "quick" run through to check that things are working?
echo "Do you wish to use a very low number of parameter configurations and epochs to quickly establish that the code is working? (y/n) "
read quick_str

if ! [[ $quick_str == "y" ||  $quick_str == "Y" || $quick_str == "n" ||  $quick_str == "N" ]]; then
    echo "Please re-run and type y or n"
    exit 1
fi

echo ""
echo "######## Starting the experiments... ############"
echo ""

source sh/Univariate.sh           # Section 2 and Section S2 of the Supplementary Material
source sh/SimulationStudies.sh    # Section 3 and Section S7 of the Supplementary Material
source sh/SimulationOnTheFly.sh   # Section S4 of the Supplementary Material
source sh/Pretraining.sh          # Section S5 of the Supplementary Material
source sh/RedSea.sh               # Section 4

echo ""
echo "######## Everything finished! ############"
echo ""
