@echo off
setlocal enabledelayedexpansion

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

for %%i in (GaussianProcess/nuFixed GaussianProcess/nuVaried Schlather ConditionalExtremes) do (
    echo.
    echo ######## Starting simulation study for %%i model ############
    echo.

    if /i "%%i" equ "GaussianProcess/nuFixed" (
        :: Generate and plot the parameter configurations
        Rscript src\%%i\Parameters.R !quick!
        :: Rscript src\%%i\PlotParameters.R # this throws an error for some reason on the HPC; it's a minor plot in the supp, so ok to omit
    )

    :: Visualise field realisations
    julia --threads=auto --project=. src\VisualiseRealisations.jl --model=%%i
    Rscript src\VisualiseRealisations.R --model=%%i

    :: Train the neural estimators
    julia --threads=auto --project=. src\Train.jl --model=%%i !quick! --m="[1]"
    julia --threads=auto --project=. src\Train.jl --model=%%i !quick! --m="[1,10,30,75,150]" --deep

    :: Plot the loss functions
    Rscript src\PlotLoss.R --path=%%i

    :: Estimation
    if /i "%%i" equ "ConditionalExtremes" (
        julia --threads=auto --project=. src\Estimate.jl --model=%%i !quick!
    ) else (
        julia --threads=auto --project=. src\Estimate.jl --model=%%i --ML !quick!
    )

    :: Results (plots, etc.)
    Rscript src\%%i\Results.R
)

endlocal
