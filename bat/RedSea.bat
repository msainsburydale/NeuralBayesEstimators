@echo off
setlocal enabledelayedexpansion

echo.
echo ######## Setting up ############
echo.

:: Check that the data is installed
if not exist data\RedSea\redseatemperature.rdata (
    echo The Red Sea data set has not been downloaded, or is in the wrong location. Please see the README for download instructions.
    exit /b 1
)

echo Do you wish to use a very low number of parameter configurations and epochs to quickly establish that the code is working? (y/n)
set /p quick_str=

if /i "!quick_str!" equ "y" (
    set quick=--quick
) elif /i "!quick_str!" equ "n" (
    set quick=
) else (
    echo Please re-run and type y or n
    exit /b 1
)

echo.
echo ######## Starting application study for the Red Sea data set ############
echo.

:: Data processing and defining the set of parameter configurations
Rscript src\RedSea\DataProcessing.R --plot
Rscript src\RedSea\Parameters.R !quick!

:: Empirical Bootstrap
julia --threads=auto --project=. src\RedSea\EmpiricalEstimatorBootstrap.jl

:: Train and estimate with each neural estimator described in the manuscript
julia --threads=auto --project=. src\RedSea\Train.jl !quick! --arch=CNN --data_type=regular
Rscript src\PlotLoss.R --path=RedSea\CNN
julia --threads=auto --project=. src\RedSea\Estimate.jl --arch=CNN --data_type=regular
Rscript src\RedSea\Results.R --arch=CNN --data_type=regular

for %%i in (regular irregular) do (
    julia --threads=auto --project=. src\RedSea\Train.jl !quick! --arch=DNN  --data_type=%%i
    Rscript src\PlotLoss.R --path=RedSea\DNN%%i
    julia --threads=auto --project=. src\RedSea\Estimate.jl --arch=DNN  --data_type=%%i
    Rscript src\RedSea\Results.R --arch=DNN --data_type=%%i
)

:: Plots that involve all of the neural estimators
Rscript src\RedSea\ResultsCombined.R

endlocal
