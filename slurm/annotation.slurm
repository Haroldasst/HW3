#!/bin/bash
#SBATCH --job-name=gtf_bed12
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/gtf_to_bed12_%j.out
#SBATCH --error=logs/gtf_to_bed12_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate snakemake

snakemake --cores "$SLURM_CPUS_PER_TASK" --sdm=conda gtf_to_bed12
