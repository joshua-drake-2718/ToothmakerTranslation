#!/usr/bin/env bash
# Build the original FORTRAN simulator (13.f90) to .tmp/fortran/toothmaker.
# Used to generate authoritative reference output for the Python translation.
set -euo pipefail

cd "$(dirname "$0")/.."

mkdir -p .tmp/fortran
gfortran -O2 -ffree-line-length-none -std=legacy \
    -J .tmp/fortran \
    -o .tmp/fortran/toothmaker \
    13.f90

echo "Built .tmp/fortran/toothmaker"
