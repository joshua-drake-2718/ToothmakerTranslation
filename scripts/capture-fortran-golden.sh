#!/usr/bin/env bash
# Run the FORTRAN simulator on examples/seal.txt and refresh the golden
# reference output in tests/golden_fortran/.
#
# Usage: scripts/capture-fortran-golden.sh
#
# Requires the binary at .tmp/fortran/toothmaker (run scripts/build-fortran.sh
# first if missing).
set -euo pipefail

cd "$(dirname "$0")/.."

if [[ ! -x .tmp/fortran/toothmaker ]]; then
    echo "Building FORTRAN binary first..."
    scripts/build-fortran.sh
fi

OUT=.tmp/fortran-baseline
rm -rf "$OUT"
mkdir -p "$OUT"

# The FORTRAN binary takes: input_file output_folder output_name iterations save_blocks
# Output filenames pad with trailing underscores due to fixed-size character buffers
# in 13.f90; the capture step renames them to clean names.
.tmp/fortran/toothmaker examples/seal.txt "$OUT" run 100 5

mkdir -p tests/golden_fortran
rm -f tests/golden_fortran/*.off

for f in "$OUT"/*.off; do
    base=$(basename "$f")
    # Strip trailing underscores in the stem, keep the .off extension.
    # e.g. "100_run_______________________.off" -> "100_run.off"
    clean=$(echo "$base" | sed -E 's/_+\.off$/.off/')
    cp "$f" "tests/golden_fortran/$clean"
done

echo "Captured FORTRAN golden output:"
ls -la tests/golden_fortran/*.off
