configfile: "config/config.yaml"

import csv
import os

with open(config["samples"]) as f:
    SAMPLES = [
        row["sample_id"]
        for row in csv.DictReader(f, delimiter="\t")
    ]

RAW_DIR = config["directories"]["raw"]
TRIMMED_DIR = config["directories"]["trimmed"]
RESULTS_DIR = config["directories"]["results"]
LOGS_DIR = config["directories"]["logs"]

INDEX_PREFIX = config["indexes"]["hisat2"]
KALLISTO_INDEX = config["indexes"]["kallisto"]
INDEX_DIR = os.path.dirname(INDEX_PREFIX)
BED12 = config["quality"]["bed12"]
RSEQC_DIR = config["quality"]["rseqc_dir"]
DEEPTOOLS_DIR = config["quality"]["deeptools_dir"]

HISAT2_INDEXES = [
    f"{INDEX_PREFIX}.{i}.ht2"
    for i in range(1, 9)
]

rule all:
    input:
        expand(
            f"{RESULTS_DIR}/fastqc/raw/{{sample}}_{{mate}}_fastqc.html",
            sample=SAMPLES, mate=[1, 2]
        ),
        expand(
            f"{RESULTS_DIR}/fastqc/raw/{{sample}}_{{mate}}_fastqc.zip",
            sample=SAMPLES, mate=[1, 2]
        ),
        expand(
            f"{RESULTS_DIR}/fastp/{{sample}}.html",
            sample=SAMPLES
        ),
        expand(
            f"{RESULTS_DIR}/fastp/{{sample}}.json",
            sample=SAMPLES
        ),
        expand(
            f"{RESULTS_DIR}/fastqc/trimmed/{{sample}}_{{mate}}_fastqc.html",
            sample=SAMPLES, mate=[1, 2]
        ),
        expand(
            f"{RESULTS_DIR}/fastqc/trimmed/{{sample}}_{{mate}}_fastqc.zip",
            sample=SAMPLES, mate=[1, 2]
        ),
        HISAT2_INDEXES,
        expand(
            f"{RESULTS_DIR}/alignment/{{sample}}.sorted.bam",
            sample=SAMPLES
        ),
        expand(
            f"{RESULTS_DIR}/alignment/{{sample}}.sorted.bam.bai",
            sample=SAMPLES
        ),
        f"{RESULTS_DIR}/featurecounts/gene_counts.txt",
        f"{RESULTS_DIR}/featurecounts/gene_counts.txt.summary",
         expand(
            f"{RESULTS_DIR}/kallisto/{{sample}}/abundance.tsv",
            sample=SAMPLES
        ),
        expand(
            f"{RESULTS_DIR}/kallisto/{{sample}}/abundance.h5",
            sample=SAMPLES
        ),
        expand(
            f"{RESULTS_DIR}/kallisto/{{sample}}/run_info.json",
            sample=SAMPLES
        ),
        BED12,
        expand(
            f"{RSEQC_DIR}/{{sample}}.infer_experiment.txt",
            sample=SAMPLES
        ),

rule fastqc_raw:
    input:
        f"{RAW_DIR}/{{sample}}_{{mate}}.fastq.gz"

    output:
        html=f"{RESULTS_DIR}/fastqc/raw/{{sample}}_{{mate}}_fastqc.html",
        zip=f"{RESULTS_DIR}/fastqc/raw/{{sample}}_{{mate}}_fastqc.zip"

    threads:
        config["threads"]["fastqc"]

    conda:
        "environment/fastqc.yaml"

    shell:
        """
        mkdir -p {RESULTS_DIR}/fastqc/raw
        fastqc --threads {threads} --outdir {RESULTS_DIR}/fastqc/raw {input}
        """


rule fastp:
    input:
        r1=f"{RAW_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{RAW_DIR}/{{sample}}_2.fastq.gz"

    output:
        r1=f"{TRIMMED_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{TRIMMED_DIR}/{{sample}}_2.fastq.gz",
        html=f"{RESULTS_DIR}/fastp/{{sample}}.html",
        json=f"{RESULTS_DIR}/fastp/{{sample}}.json"

    log:
        f"{LOGS_DIR}/fastp_{{sample}}.log"

    threads:
        config["threads"]["fastp"]

    params:
        quality=config["fastp"]["quality"],
        window_size=config["fastp"]["window_size"]

    conda:
        "environment/fastp.yaml"

    shell:
        """
        mkdir -p {TRIMMED_DIR} {RESULTS_DIR}/fastp {LOGS_DIR}

        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --detect_adapter_for_pe \
            --cut_tail \
            --cut_window_size {params.window_size} \
            --cut_mean_quality {params.quality} \
            --thread {threads} \
            --html {output.html} \
            --json {output.json} \
            > {log} 2>&1
        """


rule fastqc_trimmed:
    input:
        f"{TRIMMED_DIR}/{{sample}}_{{mate}}.fastq.gz"

    output:
        html=f"{RESULTS_DIR}/fastqc/trimmed/{{sample}}_{{mate}}_fastqc.html",
        zip=f"{RESULTS_DIR}/fastqc/trimmed/{{sample}}_{{mate}}_fastqc.zip"

    threads:
        config["threads"]["fastqc"]

    conda:
        "environment/fastqc.yaml"

    shell:
        """
        mkdir -p {RESULTS_DIR}/fastqc/trimmed
        fastqc --threads {threads} --outdir {RESULTS_DIR}/fastqc/trimmed {input}
        """

rule hisat2_index:
    input:
        genome=config["reference"]["genome"]

    output:
        HISAT2_INDEXES

    threads:
        config["threads"]["hisat2_index"]

    params:
        prefix=INDEX_PREFIX,
        index_dir=INDEX_DIR

    conda:
        "environment/alignment.yaml"

    log:
        f"{LOGS_DIR}/hisat2_index.log"

    shell:
        """
        mkdir -p {params.index_dir} {LOGS_DIR}
        hisat2-build -p {threads} {input.genome} {params.prefix} 2> {log}
        """


rule hisat2_align:
    input:
        r1=f"{TRIMMED_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{TRIMMED_DIR}/{{sample}}_2.fastq.gz",
        index=HISAT2_INDEXES

    output:
        bam=f"{RESULTS_DIR}/alignment/{{sample}}.sorted.bam"

    log:
        f"{LOGS_DIR}/hisat2/{{sample}}.log"

    threads:
        config["threads"]["hisat2_align"]

    params:
        index=INDEX_PREFIX,
        hisat2_threads=config["threads"]["hisat2_align"] - 1

    conda:
        "environment/alignment.yaml"

    shell:
        """
        mkdir -p {RESULTS_DIR}/alignment
        mkdir -p {LOGS_DIR}/hisat2

        hisat2 \
            -x {params.index} \
            -1 {input.r1} \
            -2 {input.r2} \
            -p {params.hisat2_threads} \
            2> {log} |
        samtools sort \
            -o {output.bam} -
        """

rule samtools_index:
    input:
        bam=f"{RESULTS_DIR}/alignment/{{sample}}.sorted.bam"

    output:
        bai=f"{RESULTS_DIR}/alignment/{{sample}}.sorted.bam.bai"

    threads: 1

    conda:
        "environment/alignment.yaml"

    shell:
        """
        samtools index {input.bam} {output.bai}
        """

rule featurecounts:
    input:
        bams=expand(
            f"{RESULTS_DIR}/alignment/{{sample}}.sorted.bam",
            sample=SAMPLES
        ),
        annotation=config["reference"]["annotation"]

    output:
        counts=f"{RESULTS_DIR}/featurecounts/gene_counts.txt",
        summary=f"{RESULTS_DIR}/featurecounts/gene_counts.txt.summary"

    log:
        f"{LOGS_DIR}/featurecounts.log"

    threads:
        config["threads"]["featurecounts"]

    params:
        feature_type=config["featurecounts"]["feature_type"],
        attribute=config["featurecounts"]["attribute"],
        strandedness=config["featurecounts"]["strandedness"]

    conda:
        "environment/featurecounts.yaml"

    shell:
        """
        mkdir -p {RESULTS_DIR}/featurecounts {LOGS_DIR}

        featureCounts \
            -T {threads} \
            -p \
            --countReadPairs \
            -s {params.strandedness} \
            -t {params.feature_type} \
            -g {params.attribute} \
            -a {input.annotation} \
            -o {output.counts} \
            {input.bams} \
            > {log} 2>&1
        """

rule kallisto_index:
    input:
        transcriptome=config["reference"]["transcriptome"]

    output:
        index=KALLISTO_INDEX

    threads:
        config["threads"]["kallisto_index"]

    params:
        index_dir=os.path.dirname(KALLISTO_INDEX)

    conda:
        "environment/kallisto.yaml"

    log:
        f"{LOGS_DIR}/kallisto_index.log"

    shell:
        """
        mkdir -p {params.index_dir} {LOGS_DIR}

        kallisto index \
            -i {output.index} \
            -t {threads} \
            {input.transcriptome} \
            > {log} 2>&1
        """

rule kallisto_quant:
    input:
        index=KALLISTO_INDEX,
        r1=f"{TRIMMED_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{TRIMMED_DIR}/{{sample}}_2.fastq.gz"

    output:
        abundance=f"{RESULTS_DIR}/kallisto/{{sample}}/abundance.tsv",
        h5=f"{RESULTS_DIR}/kallisto/{{sample}}/abundance.h5",
        run_info=f"{RESULTS_DIR}/kallisto/{{sample}}/run_info.json"

    log:
        f"{LOGS_DIR}/kallisto/{{sample}}.log"

    threads:
        config["threads"]["kallisto_quant"]

    params:
        outdir=f"{RESULTS_DIR}/kallisto/{{sample}}"

    conda:
        "environment/kallisto.yaml"

    shell:
        """
        mkdir -p {params.outdir} {LOGS_DIR}/kallisto

        kallisto quant \
            -i {input.index} \
            -o {params.outdir} \
            -t {threads} \
            {input.r1} \
            {input.r2} \
            > {log} 2>&1
        """

rule gtf_to_bed12:
    input:
        gtf=config["reference"]["annotation"]

    output:
        genepred="reference/index/gencode.v50.annotation.genePred",
        bed12=BED12

    log:
        f"{LOGS_DIR}/gtf_to_bed12.log"

    conda:
        "environment/rseqc.yaml"

    shell:
        """
        mkdir -p reference/index {LOGS_DIR}

        gtfToGenePred {input.gtf} {output.genepred} 2> {log}

        genePredToBed {output.genepred} {output.bed12} >> {log} 2>&1
        """
rule infer_experiment:
    input:
        bam=f"{RESULTS_DIR}/alignment/{{sample}}.sorted.bam",
        bed12=BED12

    output:
        report=f"{RSEQC_DIR}/{{sample}}.infer_experiment.txt"

    log:
        f"{LOGS_DIR}/infer_experiment/{{sample}}.log"

    conda:
        "environment/rseqc.yaml"

    shell:
        """
        mkdir -p {RSEQC_DIR} $(dirname {log})

        infer_experiment.py \
            -r {input.bed12} \
            -i {input.bam} \
            > {output.report} 2> {log}
        """
