#!/bin/bash

# Run async.py on everything but SILHS code
find ../src -maxdepth 1 -path "./src/SILHS" -prune -o -name "*.F90" -print0 \
  | xargs -P 12 -0 -I {} python async.py {}