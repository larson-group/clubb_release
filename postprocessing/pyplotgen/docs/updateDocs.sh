#!/usr/bin/env bash

echo "Updating HTML Reference Documentation for PyPlotgen"
echo "-------------------------"

# if there's a venv, try to use it
source ../venv/bin/activate

echo "Ensuring python dependencies are installed"
pip3 install sphinx sphinx_rtd_theme

make clean
make html

echo "If the build succeeded, open '_build/html/index.html' to view the reference docs."
