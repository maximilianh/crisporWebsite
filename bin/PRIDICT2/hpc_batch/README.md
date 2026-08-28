# Running PRIDICT2.0 in Batch Mode on HPC with SLURM

This guide outlines the process of running PRIDICT2.0 in batch mode on a High-Performance Computing (HPC) cluster using the SLURM job management system.

## Prerequisites

- Access to an HPC cluster with SLURM job management
- PRIDICT2.0 environment installed on the cluster (refer to the [main PRIDICT2.0 installation guide](../README.md#installation-using-anaconda-linux-or-mac-os-) for instructions)
- Basic knowledge of command-line interface and SLURM commands

## Steps

Move into the hpc_batch folder of the repository:
```bash
cd hpc_batch
```

### 1. Prepare Batch Files

Two shell scripts are provided for batch processing:

- `cluster_part1_split_files.sh`: Splits the input file into smaller parts
- `cluster_part2_slurm_job.sh`: Submits SLURM jobs for each part

Make these scripts executable:

```bash
chmod +x cluster_part1_split_files.sh cluster_part2_slurm_job.sh
```

### 2. Split Input File

Decide how many splits you want to make (num_splits) and run:

```bash
./cluster_part1_split_files.sh --input_folder "../input/" --num_splits <number_of_splits> <input_file>
```

Replace `<number_of_splits>` with your desired number and `<input_file>` with your input file name (e.g. `batch_template.csv`).

### 3. Modify SLURM Job Script

Edit the `cluster_part2_slurm_job.sh` file to adjust the `--array` parameter to match the number of splits:

1. Open the file in a text editor:
   ```bash
   nano cluster_part2_slurm_job.sh
   ```

2. Find the line starting with `#SBATCH --array=` and modify it to:
   ```bash
   #SBATCH --array=0-<num_splits-1>
   ```
   Replace `<num_splits-1>` with the number of splits you created minus one.

3. Save and exit the editor.

### 4. Submit SLURM Jobs

After splitting the input file and modifying the job script, submit the SLURM jobs:

```bash
sbatch ./cluster_part2_slurm_job.sh ~/data/PRIDICT2 batch_template
```

Adjust the path and base name (batch filename without extension) as necessary.

### 5. Monitor Jobs

Use SLURM commands to monitor your jobs:

```bash
squeue -u your_username
```

### 6. Collect Results

Once all jobs are complete, collect and combine the results from the `predictions` directory.

## Troubleshooting

- If jobs fail, check the output files (`job_pridict2_*.out`) for error messages.
- Ensure you have sufficient permissions and resources allocated for your jobs.
- Verify that the PRIDICT2.0 environment is correctly set up and accessible to the compute nodes.

## Additional Notes

- Adjust other SLURM parameters in `cluster_part2_slurm_job.sh` (e.g., `--cpus-per-task`, `--mem-per-cpu`, `--time`) based on your cluster's policies and the requirements of your specific batch job.

For further assistance or questions, please contact your cluster's support team or open an issue on the PRIDICT2.0 GitHub repository.