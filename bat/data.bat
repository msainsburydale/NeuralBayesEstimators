@echo off
setlocal enabledelayedexpansion

if exist data\RedSea\redseatemperature.rdata (
    echo The Red Sea data has already been downloaded.
) else (
    rem Use curl to download the file (assuming curl is available in your environment)
    curl -o data\RedSea\redseatemperature.rdata https://zenodo.org/record/8134200/files/redseatemperature.rdata
)
