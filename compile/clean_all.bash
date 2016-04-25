#!/bin/bash
# $Id: compile.bash 5079 2011-04-01 17:25:38Z dschanen@uwm.edu $

# What directory am I located in?
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd "$DIR/../bin" && make distclean
