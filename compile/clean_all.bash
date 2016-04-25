#!/bin/bash
# $Id$

# What directory am I located in?
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd "$DIR/../bin" && make distclean
