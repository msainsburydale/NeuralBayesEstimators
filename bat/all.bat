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

:: Should we do a "quick" run through to check that things are working?
set /p quick_str=Do you wish to use a very low number of parameter configurations and epochs to quickly establish that the code is working? (y/n)
if /i not "!quick_str!"=="y" if /i not "!quick_str!"=="n" (
    echo Please re-run and type y or n
    exit /b 1
)

echo.
echo ######## Starting the experiments... ############
echo.

call sh\Univariate.bat          :: Section 2 and Section S2 of the Supplementary Material
call sh\SimulationStudies.bat   :: Section 3 and Section S7 of the Supplementary Material
call sh\SimulationOnTheFly.bat  :: Section S4 of the Supplementary Material
call sh\Pretraining.bat         :: Section S5 of the Supplementary Material
call sh\RedSea.bat              :: Section 4

echo.
echo ######## Everything finished! ############
echo.
