#!/bin/bash

# Default input folder
INPUT_FOLDER="../input/"

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input_folder) INPUT_FOLDER="$2"; shift ;;
        --num_splits) NUM_SPLITS="$2"; shift ;;
        *) INPUT_FILE="$1" ;;
    esac
    shift
done

# Check if input file was provided
if [ -z "$INPUT_FILE" ]; then
    echo "Error: No input file provided."
    echo "Usage: $0 [--input_folder <folder>] --num_splits <number> <input_file>"
    exit 1
fi

# Check if num_splits was provided
if [ -z "$NUM_SPLITS" ]; then
    echo "Error: --num_splits argument is required."
    echo "Usage: $0 [--input_folder <folder>] --num_splits <number> <input_file>"
    exit 1
fi

# Construct full input path
FULL_INPUT_PATH="${INPUT_FOLDER%/}/${INPUT_FILE}"

# Check if file exists
if [ ! -f "$FULL_INPUT_PATH" ]; then
    echo "Error: File '$FULL_INPUT_PATH' not found."
    exit 1
fi

# Variables
HEADER=$(head -1 "$FULL_INPUT_PATH")

# Calculate the number of lines per split file
TOTAL_LINES=$(($(wc -l < "$FULL_INPUT_PATH") - 1)) # Exclude header
LINES_PER_FILE=$((TOTAL_LINES / NUM_SPLITS))

# Adjust if LINES_PER_FILE is less than 1
if [ "$LINES_PER_FILE" -lt 1 ]; then
    LINES_PER_FILE=1
fi

# Split the file, preserving the header in each part
tail -n +2 "$FULL_INPUT_PATH" | split -l $LINES_PER_FILE -d - "${FULL_INPUT_PATH%.*}_part"

# Add the header to each split file and rename them
COUNTER=0
for FILE in ${FULL_INPUT_PATH%.*}_part*; do
    NEW_NAME="${FULL_INPUT_PATH%.*}_${COUNTER}.csv"
    mv "$FILE" "$NEW_NAME"
    sed -i "1s/^/${HEADER}\n/" "$NEW_NAME"
    COUNTER=$((COUNTER + 1))
done
