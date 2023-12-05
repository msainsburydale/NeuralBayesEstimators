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

for %%model in (GaussianProcess\nuFixed GaussianProcess\nuVaried Schlather ConditionalExtremes) do (
    echo.
    echo ######## Starting simulation study for %%model model ############
    echo.

    if /i "%%model"=="GaussianProcess\nuFixed" (
        :: Generate and plot the parameter configurations
        Rscript src\%%model\Parameters.R !quick!
        :: Rscript src\%%model\PlotParameters.R # this throws an error for some reason on the HPC; it's a minor plot in the supp, so ok to omit
    )

    :: Visualise field realisations
    julia --threads=auto --project=. src\VisualiseRealisations.jl --model=%%model
    Rscript src\VisualiseRealisations.R --model=%%model

    :: Train the neural estimators
    julia --threads=auto --project=. src\Train.jl --model=%%model !quick! --m="[1]"
    julia --threads=auto --project=. src\Train.jl --model=%%model !quick! --m="[1,10,30,75,150]" --deep

    :: Plot the loss functions
    Rscript src\PlotLoss.R --path=%%model

    :: Estimation
    if /i "%%model"=="ConditionalExtremes" (
        julia --threads=auto --project=. src\Estimate.jl --model=%%model !quick!
    ) else (
        julia --threads=auto --project=. src\Estimate.jl --model=%%model --ML !quick!
    )

    :: Results (plots, etc.)
    Rscript src\%%model\Results.R
)
