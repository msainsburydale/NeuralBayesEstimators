#!/bin/bash

set -e

echo "Really clear all intermediates and results? (y/n)"
read answer

if [[ $answer == "y" || $answer == "Y" ]]; then
    rm -rf intermediates/*
    rm -rf img/*
    rm -rf results/*
fi




