#!/bin/bash

# Check if at least one argument is provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 (<fib-file> | <nodes> <connections> <output-filename>) [c_thr] [start_i] [numcycles] [bell] [direction]"
    exit 1
fi

# Assign parameters with defaults
C_THR="${4:-0.9}"
START_I="${5:-10}"
NUMCYCLES="${6:-10}"
BELL="${7:-5}"
DIRECTED="${8:-1}"

# Define paths
BUNDLER="./bundler/build/Desktop-Debug/bundler"
FIBVIEWER="./fibviewer/build/Desktop-Debug/fibviewer"
OUTPUT_DIR="./outputs"

# Ensure binaries are executable
if [[ ! -x "$BUNDLER" ]]; then
    chmod +x "$BUNDLER"
fi

if [[ ! -x "$FIBVIEWER" ]]; then
    chmod +x "$FIBVIEWER"
fi

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Determine whether to use -fib OR -nodes & -cons
if [[ "$1" == *.fib ]]; then
    # Case: Using a .fib file
    FIB_FILE="$1"
    echo "Running Bundler with: -fib $FIB_FILE -c_thr $C_THR -start_i $START_I -numcycles $NUMCYCLES"
    "$BUNDLER" -fib "$FIB_FILE" -c_thr "$C_THR" -start_i "$START_I" -numcycles "$NUMCYCLES" -directed "$DIRECTED"

    # Move and run FibViewer on the output file
    FIB_TXT_FILE="${FIB_FILE}_c_thr$(printf "%.4f" "$C_THR")_numcycles$(printf "%02d" "$NUMCYCLES")_start_i$(printf "%04d" "$START_I").vtk"
else
    # Case: Using -nodes and -cons
    if [ "$#" -lt 3 ]; then
        echo "Error: If using -nodes, you must provide <nodes> <connections> and <output-filename>!"
        exit 1
    fi

    NODES="$1"
    CONNECTIONS="$2"
    OUTPUT_FILENAME="$3"

    echo "Running Bundler with: -nodes $NODES -cons $CONNECTIONS -fileName $OUTPUT_FILENAME -c_thr $C_THR -start_i $START_I -numcycles $NUMCYCLES"
    "$BUNDLER" -nodes "$NODES" -cons "$CONNECTIONS" -fileName "$OUTPUT_FILENAME" -c_thr "$C_THR" -start_i "$START_I" -numcycles "$NUMCYCLES" -directed "$DIRECTED"

    # Move and run FibViewer on the output file
    FIB_TXT_FILE="${OUTPUT_FILENAME}_c_thr$(printf "%.4f" "$C_THR")_numcycles$(printf "%02d" "$NUMCYCLES")_start_i$(printf "%04d" "$START_I").vtk"
fi

# Check if the output file was created
if [[ -f "$FIB_TXT_FILE" ]]; then
    # Move output file to the outputs folder
    mv "$FIB_TXT_FILE" "$OUTPUT_DIR/"
    echo "Moved output file to $OUTPUT_DIR/$FIB_TXT_FILE"

    # Run FibViewer on the moved file
    echo "Running FibViewer with: $OUTPUT_DIR/$FIB_TXT_FILE"

    python3 ./metrics/voxel_based_edge_density.py "$OUTPUT_DIR/$FIB_TXT_FILE"

    "$FIBVIEWER" "$OUTPUT_DIR/$FIB_TXT_FILE"
else
    echo "Error: Expected output file '$FIB_TXT_FILE' not found!"
    exit 2
fi
