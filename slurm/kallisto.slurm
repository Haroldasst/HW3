#!/bin/bash

#SBATCH --job-name=HW3_Kallisto
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/kallisto_%j.out
#SBATCH --error=logs/kallisto_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate snakemake

snakemake --cores "$SLURM_CPUS_PER_TASK" --sdm=conda
