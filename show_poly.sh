#!/bin/bash

FIBVIEWER="./fibviewer/build/Desktop-Debug/fibviewer"
VIEWFILE="./viewmatrix.txt"

# Loop through files that start with "poly_test"
for test_file in test*; do
    # Check if it's actually a file
    if [[ -f "$test_file" ]]; then
        python3 ./metrics/voxel_based_edge_density.py "$test_file"
        "$FIBVIEWER" "$test_file" -screenshot &
    fi
done

wait
