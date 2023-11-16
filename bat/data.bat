@echo off
setlocal enabledelayedexpansion

if exist data\RedSea\redseatemperature.rdata (
    echo The Red Sea data has already been downloaded.
) else (
    :: Use curl for Windows to download the file
    curl -o data\RedSea\redseatemperature.rdata https://zenodo.org/record/8134200/files/redseatemperature.rdata
    :: Alternatively, you can use wget if it's available on your system
    :: wget -P data/RedSea https://zenodo.org/record/8134200/files/redseatemperature.rdata
)

endlocal
