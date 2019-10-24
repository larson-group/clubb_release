#!/usr/bin/env bash

echo "Updating HTML Reference Documentation for pyplotgen"
echo "-------------------------"

# if there's a venv, try to use it
source ../venv/bin/activate

sphinx-build -b html ./config ./html_docs -a

