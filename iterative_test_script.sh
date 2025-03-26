#!/bin/bash

# Check if at least one argument is provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 (<fib-file> | <nodes> <connections> <output-filename>) [c_thr] [start_i] [numcycles] [bell]"
    exit 1
fi

# Assign parameters with defaults (upper bounds for start_i and numcycles)
C_THR="${4:-0.9}"
MAX_START_I="${5:-10}"     # Upper bound for start_i
MAX_NUMCYCLES="${6:-10}"   # Upper bound for numcycles
BELL="${7:-5}"

# Define paths
BUNDLER="./bundler/build/Desktop-Debug/bundler"
FIBVIEWER="./fibviewer/build/Desktop-Debug/fibviewer"
OUTPUT_DIR="./outputs"

# Ensure binaries are executable
chmod +x "$BUNDLER" 2>/dev/null
chmod +x "$FIBVIEWER" 2>/dev/null

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Initialize counters
current_start_i=1
current_numcycles=1

# Function to run Bundler and FibViewer
run_processing() {
    local start_i=$1
    local numcycles=$2

    if [[ "$FIB_FILE" ]]; then
        echo "Running Bundler with: -fib $FIB_FILE -c_thr $C_THR -start_i $start_i -numcycles $numcycles"
        "$BUNDLER" -fib "$FIB_FILE" -c_thr "$C_THR" -start_i "$start_i" -numcycles "$numcycles" &
        local output_file="${FIB_FILE}.vtk"
    else
        echo "Running Bundler with: -nodes $NODES -cons $CONNECTIONS -fileName $OUTPUT_FILENAME -c_thr $C_THR -start_i $start_i -numcycles $numcycles"
        "$BUNDLER" -nodes "$NODES" -cons "$CONNECTIONS" -fileName "$OUTPUT_FILENAME" -c_thr "$C_THR" -start_i "$start_i" -numcycles "$numcycles" &
        local output_file="${OUTPUT_FILENAME}.vtk"
    fi

    # Wait for Bundler to finish
    wait

    # Check if output file exists
    if [[ -f "$output_file" ]]; then
        mv "$output_file" "$OUTPUT_DIR/"
        echo "Moved output file to $OUTPUT_DIR/$output_file"

        # Run FibViewer in the background
        echo "Running FibViewer with: $OUTPUT_DIR/$output_file"
        "$FIBVIEWER" "$OUTPUT_DIR/$output_file" &
    else
        echo "Error: Expected output file '$output_file' not found!"
    fi
}

# Determine whether to use -fib OR -nodes & -cons
if [[ "$1" == *.fib ]]; then
    FIB_FILE="$1"
else
    if [ "$#" -lt 3 ]; then
        echo "Error: If using -nodes, you must provide <nodes> <connections> and <output-filename>!"
        exit 1
    fi
    NODES="$1"
    CONNECTIONS="$2"
    OUTPUT_FILENAME="$3"
fi

# Loop to execute in parallel with cycling start_i and increasing numcycles
while [[ $current_numcycles -le $MAX_NUMCYCLES ]]; do
    run_processing "$current_start_i" "$current_numcycles"

    ((current_start_i++))
    if [[ $current_start_i -gt $MAX_START_I ]]; then
        current_start_i=1
        ((current_numcycles++))
    fi

    # Ensure we do not exceed max numcycles
    if [[ $current_numcycles -gt $MAX_NUMCYCLES ]]; then
        break
    fi
done

# Wait for all background processes to complete
wait
echo "Processing complete."
