#!/bin/bash

# Paths
BUNDLER="./bundler/build/Desktop-Debug/bundler"
FIBVIEWER="./fibviewer/build/Desktop-Debug/fibviewer"
BUNDLER_OUTPUT_DIR="./"
SCREENSHOT_DIR="./outputs/screenshots"
TMP_VTK_DIR="./outputs/tmp_vtks"
FINAL_OUTPUT_DIR="./outputs/final"
PYTHON_SCRIPT="./make_gif_or_slider.py"
VIEWFILE="./viewmatrix.txt"

# Clean previous screenshots
rm -rf "$SCREENSHOT_DIR"
rm -rf "$TMP_VTK_DIR"
mkdir -p "$SCREENSHOT_DIR"
mkdir -p "$FINAL_OUTPUT_DIR"
mkdir -p "$TMP_VTK_DIR"

# Run bundler
echo "Running bundler..."

# Check if at least one argument is provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 (<fib-file> | <nodes> <connections> <output-filename>) [c_thr] [start_i] [numcycles] [bell]"
    exit 1
fi

# Assign parameters with defaults
C_THR="${4:-0.9}"
START_I="${5:-10}"
NUMCYCLES="${6:-10}"
BELL="${7:-5}"

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
    "$BUNDLER" -fib "$FIB_FILE" -c_thr "$C_THR" -start_i "$START_I" -numcycles "$NUMCYCLES" -checkpoints 1 -directed 1

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
    "$BUNDLER" -nodes "$NODES" -cons "$CONNECTIONS" -fileName "$OUTPUT_FILENAME" -c_thr "$C_THR" -start_i "$START_I" -numcycles "$NUMCYCLES" -checkpoints 1 -directed 1

    # Move and run FibViewer on the output file
    FIB_TXT_FILE="${OUTPUT_FILENAME}_c_thr$(printf "%.4f" "$C_THR")_numcycles$(printf "%02d" "$NUMCYCLES")_start_i$(printf "%04d" "$START_I").vtk"
fi


# Find all VTK files generated that match the OUTPUT_FILENAME prefix
vtk_files=$(find "$BUNDLER_OUTPUT_DIR" -name "$OUTPUT_FILENAME*.vtk" | sort)

screenshot_index=0

# Loop through each .vtk file
for vtk_file in $vtk_files; do
    echo "Processing $vtk_file..."

    # Ensure that the file exists and is being moved correctly
    if [[ ! -f "$vtk_file" ]]; then
        echo "Error: VTK file does not exist: $vtk_file"
        continue
    fi

    # Move the .vtk file to the dedicated tmp_vtks folder
    mv "$vtk_file" "$TMP_VTK_DIR/"

    # Extract the basename of the vtk file for later use
    vtk_basename=$(basename "$vtk_file")

    # Run FibViewer, it will generate a screenshot with the same name (vtk_basename.vtk.png)
    echo "$FIBVIEWER $TMP_VTK_DIR/$vtk_basename -screenshot -viewmatrix $VIEWFILE"
    "$FIBVIEWER" "$TMP_VTK_DIR/$vtk_basename" -screenshot -viewmatrix "$VIEWFILE"

    # The screenshot will be generated with the same name as the .vtk file but with .vtk.png appended
    screenshot_file="$TMP_VTK_DIR/$vtk_basename.png"

    # After FibViewer is done, move the screenshot to the main screenshots directory
    mv "$screenshot_file" "$SCREENSHOT_DIR"

    ((screenshot_index++))
done

# Wait if needed
wait

# Create GIF or slideshow
echo "Creating animation..."
python3 "$PYTHON_SCRIPT" "$SCREENSHOT_DIR" "$FINAL_OUTPUT_DIR"

echo "All done!"
