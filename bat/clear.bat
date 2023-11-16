@echo off
setlocal enabledelayedexpansion

echo Really clear all intermediates and results? (y/n)
set /p answer=

if /i "!answer!" equ "y" (
    rmdir /s /q intermediates
    rmdir /s /q img
    rmdir /s /q results
)

endlocal
