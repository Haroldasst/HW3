HW3 - Eukaryotic RNA-seq analysis using Snakemake

1. PROJECT DESCRIPTION

Aim - learn to develop reproducible bioinformatics workflows
using Snakemake, perform alignment-based and alignment-free 
RNA-seq quantification, evaluate sequencing and alignment quality,
and use Conda, Git, and SLURM in a bioinformatics project.

2. DATASET

GEO accession: GSE52778

Samples:
GSM1275862 - Untreated
GSM1275863 - Dexamethasone
GSM1275864 - Albuterol
GSM1275866 - Untreated

SRA accession numbers:
SRR1039508
SRR1039509
SRR1039510
SRR1039512

Organism: Homo sapiens
Reference genome: GRCh38

3. SOFTWARE

- Snakemake
- Conda
- FastQC
- fastp
- HISAT2
- SAMtools
- featureCounts
- kallisto
- RSeQC
- deepTools
- MultiQC
- IGV

4. PROJECT STRUCTURE

config/ - configuration files and sample information
envs/ - Conda environment definitions
scripts/ - additional scripts
data/ - sequencing data
reference/ - reference genome and annotation
results/ - analysis results
logs/ - program logs
slurm/ - SLURM job scripts

5. WORKFLOW EXECUTION

- Downloading sequencing data.
- Raw read quality control using FastQC.
- Adapter and quality trimming using fastp.
- HISAT2 alignment and featureCounts quantification.
- Transcript quantification using kallisto.
- Alignment quality assessment using RSeQC.
- Sample correlation and PCA using deepTools.
- Quality control reporting using MultiQC.

6. WORKFLOW EXECUTION

Check the workflow:

snakemake -n --cores 4 --sdm=conda

Run the workflow on HPC using SLURM:

sbatch slurm/<script_name>.slurm

Sample information and analysis parameters can be changed in the config/ directory.

7. RESULTS

Analysis results are stored in the results/ directory.

The project also includes FastQC reports, featureCounts results, kallisto outputs, RSeQC and deepTools plots, and a MultiQC report.
