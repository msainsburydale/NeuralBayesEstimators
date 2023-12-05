@echo off
setlocal enabledelayedexpansion

if not defined quick_str (
    echo Do you wish to use a very low number of parameter configurations and epochs to quickly establish that the code is working? (y/n)
    set /p quick_str=
)

if /i "!quick_str!"=="y" (
    set "quick=--quick"
    set "ML="
) else if /i "!quick_str!"=="n" (
    set "quick="
    set "ML=--ML"
) else (
    echo Please re-run and type y or n
    exit /b 1
)

echo.
echo ######## Starting experiment on the effect of pretraining ############
echo.

julia --threads=auto --project=. src\Pretraining\Train.jl !quick!
Rscript src\Pretraining\Results.R
