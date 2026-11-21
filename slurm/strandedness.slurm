#!/bin/bash
#SBATCH --job-name=strandedness
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/strandedness_%j.out
#SBATCH --error=logs/strandedness_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate snakemake

snakemake --cores "$SLURM_CPUS_PER_TASK" --sdm=conda
