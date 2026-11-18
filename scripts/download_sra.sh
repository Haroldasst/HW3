#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

SAMPLES="config/samples.tsv"
SRA_DIR="data/sra"
RAW_DIR="data/raw"
TMP_DIR="data/tmp"

THREADS="$SLURM_CPUS_PER_TASK"

mkdir -p "$SRA_DIR" "$RAW_DIR" "$TMP_DIR"

tail -n +2 "$SAMPLES" |
while IFS=$'\t' read -r sample_id geo_accession sra_accession cell_line treatment
do
    echo "Processing $sample_id"

    if [[ -s "$RAW_DIR/${sample_id}_1.fastq.gz" && 
		  -s "$RAW_DIR/${sample_id}_2.fastq.gz" ]]; then
        echo "Already exists: $sample_id"
        continue
    fi

    prefetch --output-directory "$SRA_DIR" "$sra_accession"

    SRA_FILE="$SRA_DIR/$sra_accession/$sra_accession.sra"

    vdb-validate "$SRA_FILE"

    mkdir -p "$TMP_DIR/$sra_accession"

    fasterq-dump "$SRA_FILE" \
        --split-files \
        --threads "$THREADS" \
        --outdir "$RAW_DIR" \
        --temp "$TMP_DIR/$sra_accession"

    R1="$RAW_DIR/${sra_accession}_1.fastq"
    R2="$RAW_DIR/${sra_accession}_2.fastq"

    if [[ ! -s "$R1" || ! -s "$R2" ]]; then
        echo "ERROR: FASTQ files missing for $sample_id"
        exit 1
    fi

    gzip -f "$R1" "$R2"

    mv "$R1.gz" "$RAW_DIR/${sample_id}_1.fastq.gz"
    mv "$R2.gz" "$RAW_DIR/${sample_id}_2.fastq.gz"

    rm -rf "$TMP_DIR/$sra_accession"

    echo "Completed $sample_id"

done
