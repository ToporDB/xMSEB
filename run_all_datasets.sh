#!/bin/bash

INPUT_DIR="./inputs"
SCRIPT="./test_script.sh"
C_THR="0.8"
START_I="40"
NUMCYCLES="10"
BELL="5"
DIRECTED="0"

# Loop through directories in the inputs folder
for dataset_dir in "$INPUT_DIR"/*/; do
    # Check for nodes and edges files
    nodes_file=$(find "$dataset_dir" -type f -name "*nodes.txt" | head -n 1)
    edges_file=$(find "$dataset_dir" -type f -name "*edges.txt" | head -n 1)

    if [[ -f "$nodes_file" && -f "$edges_file" ]]; then
        # Extract a clean dataset name
        dataset_name=$(basename "$dataset_dir" | sed 's/[^a-zA-Z0-9]/_/g')

        echo "Processing dataset: $dataset_name"
        echo "Nodes file: $nodes_file"
        echo "Edges file: $edges_file"

        # Run the bundling script
        "$SCRIPT" "$nodes_file" "$edges_file" "$dataset_name" "$C_THR" "$START_I" "$NUMCYCLES" "$BELL" "$DIRECTED"
    else
        echo "Skipping $dataset_dir — nodes or edges file not found."
    fi
done

echo "All datasets processed."

