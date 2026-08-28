#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=12:00:00
#SBATCH --output=job_pridict2_%a.out
#SBATCH --array=0-10

# Check if the required arguments are provided
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Error: Missing arguments."
    echo "Usage: sbatch $0 <root_folder> <base_name>"
    exit 1
fi

# Assign the arguments to variables
ROOT_FOLDER=$1
BASE_NAME=$2

# Change to the root directory
cd "$ROOT_FOLDER" || { echo "Error: Failed to change to directory $ROOT_FOLDER"; exit 1; }

# Get the file part number from the array job
FILE_PART=$(($SLURM_ARRAY_TASK_ID))

# Run the Python script with the specific part of the CSV
conda run -n pridict2 python3 "${ROOT_FOLDER}/pridict2_pegRNA_design.py" batch \
    --input-dir "${ROOT_FOLDER}/input" \
    --input-fname "${BASE_NAME}_${FILE_PART}.csv" \
    --output-dir "${ROOT_FOLDER}/predictions" \
    --cores 1