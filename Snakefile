import csv

with open("config/samples.tsv") as f:
    SAMPLES = [
        row["sample_id"]
        for row in csv.DictReader(f, delimiter="\t")
    ]


rule all:
    input:
        expand(
            "results/fastqc/raw/{sample}_{mate}_fastqc.html",
            sample=SAMPLES,
            mate=[1, 2]
        ),
        expand(
            "results/fastqc/raw/{sample}_{mate}_fastqc.zip",
            sample=SAMPLES,
            mate=[1, 2]
        )


rule fastqc_raw:
    input:
        "data/raw/{sample}_{mate}.fastq.gz"

    output:
        html="results/fastqc/raw/{sample}_{mate}_fastqc.html",
        zip="results/fastqc/raw/{sample}_{mate}_fastqc.zip"

    threads: 2

    conda:
        "environment/fastqc.yaml"

    shell:
        """
        mkdir -p results/fastqc/raw
        fastqc --threads {threads} --outdir results/fastqc/raw {input}
        """
