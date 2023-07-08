#!/bin/bash

set -e

if [ -e  data/RedSea/redseatemperature.rdata ]
then
  echo "The Red Sea data has already been downloaded."
else
  wget --directory=data/RedSea https://hpc.niasra.uow.edu.au/ckan/dataset/ccbae377-d180-4441-b27e-c015bf5f89e2/resource/1484d341-9b0d-4b4a-8b4a-4b61ea4d236f/download/redseatemperature.rdata
fi
