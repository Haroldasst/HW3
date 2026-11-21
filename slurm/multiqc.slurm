#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/multiqc_%j.out
#SBATCH --error=logs/multiqc_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate snakemake

snakemake --cores "$SLURM_CPUS_PER_TASK" --sdm=conda
