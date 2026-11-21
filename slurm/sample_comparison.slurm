#!/bin/bash
#SBATCH --job-name=sample_comparison
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/sample_comparison_%j.out
#SBATCH --error=logs/sample_comparison_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate snakemake

snakemake --cores "$SLURM_CPUS_PER_TASK" --sdm=conda
