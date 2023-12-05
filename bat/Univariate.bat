@echo off
setlocal enabledelayedexpansion

if not defined quick_str (
    echo Do you wish to use a very low number of parameter configurations and epochs to quickly establish that the code is working? (y/n)
    set /p quick_str=
)

if /i "!quick_str!"=="y" (
    set "quick=--quick"
) else if /i "!quick_str!"=="n" (
    set "quick="
) else (
    echo Please re-run and type y or n
    exit /b 1
)

set "model=Uniform"

echo.
echo ######## Starting experiment on approximating the Bayes estimator: %model% distribution ############
echo.

:: Show that NN estimators can approximate the Bayes estimator for a one-parameter model
julia --threads=auto --project=. src\Univariate\Train.jl --model=%model% !quick! --aggregation=logsumexp --m=10
julia --threads=auto --project=. src\Univariate\Estimate.jl --model=%model%
Rscript src\Univariate\ResultsFixedSampleSize.R --model=%model%

set "model=Normalsigma"

echo.
echo ######## Starting experiment on the effect of variable sample size: %model% distribution ############
echo.

:: Train the estimators
julia --threads=auto --project=. src\Univariate\Train.jl --model=%model% !quick! --m=5
julia --threads=auto --project=. src\Univariate\Train.jl --model=%model% !quick! --m=150
julia --threads=auto --project=. src\Univariate\Train.jl --model=%model% !quick! --m=-150

:: Estimate and plot results
julia --threads=auto --project=. src\Univariate\Estimate.jl --model=%model%
Rscript src\Univariate\ResultsVariableSampleSize.R --model=%model%
