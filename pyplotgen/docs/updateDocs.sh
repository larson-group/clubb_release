#!/usr/bin/env bash

echo "Updating HTML Reference Documentation for PyPlotgen"
echo "-------------------------"

pip install --user virtualenv
virtualenv sphinx-venv

# if there's a venv, try to use it
source sphinx-venv/bin/activate

echo "Ensuring python dependencies are installed"
pip3 install -r ../requirements.txt
pip3 install sphinx sphinx_rtd_theme

mkdir _static
mkdir _templates

sphinx-build -o . ..

make clean
make html

echo "If the build succeeded, open '_build/html/index.html' to view the reference docs."
deactivate