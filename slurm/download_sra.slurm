#!/bin/bash

#SBATCH --job-name=HW3_SRA
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/download_%j.out
#SBATCH --error=logs/download_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

bash scripts/download_sra.sh

