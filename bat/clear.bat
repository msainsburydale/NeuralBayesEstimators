@echo off
setlocal enabledelayedexpansion

echo Really clear all intermediates and results? (y/n)
set /p answer=

if /i "!answer!"=="y" (
    del /q intermediates\* 2>nul
    del /q img\* 2>nul
    del /q results\* 2>nul
)
