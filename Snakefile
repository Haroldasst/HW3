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
BIGWIG_DIR = f"{DEEPTOOLS_DIR}/bigwig"
GENEBODY_DIR = f"{RSEQC_DIR}/gene_body_coverage"
BED12_SUBSET = config["quality"]["bed12_subset"]
INNER_DISTANCE_DIR = f"{RSEQC_DIR}/inner_distance"
COMPARISON_DIR = f"{DEEPTOOLS_DIR}/sample_comparison"

BIGWIG_SUMMARY = f"{COMPARISON_DIR}/multiBigwigSummary.npz"
BIGWIG_RAW = f"{COMPARISON_DIR}/multiBigwigSummary.tab"

CORRELATION_PDF = f"{COMPARISON_DIR}/correlation_heatmap.pdf"
CORRELATION_MATRIX = f"{COMPARISON_DIR}/correlation_matrix.tsv"

PCA_PDF = f"{COMPARISON_DIR}/pca.pdf"
PCA_DATA = f"{COMPARISON_DIR}/pca_data.tsv"

SAMPLE_COMPARISON_OUTPUTS = [
    BIGWIG_SUMMARY,
    BIGWIG_RAW,
    CORRELATION_PDF,
    CORRELATION_MATRIX,
    PCA_PDF,
    PCA_DATA
]

INNER_DISTANCE_OUTPUTS = expand(
    [
        f"{INNER_DISTANCE_DIR}/{{sample}}.inner_distance.txt",
        f"{INNER_DISTANCE_DIR}/{{sample}}.inner_distance_freq.txt",
        f"{INNER_DISTANCE_DIR}/{{sample}}.inner_distance_plot.pdf",
        f"{INNER_DISTANCE_DIR}/{{sample}}.inner_distance_plot.r"
    ],
    sample=SAMPLES
)

GENEBODY_OUTPUTS = expand(
    [
        f"{GENEBODY_DIR}/{{sample}}.geneBodyCoverage.txt",
        f"{GENEBODY_DIR}/{{sample}}.geneBodyCoverage.pdf",
        f"{GENEBODY_DIR}/{{sample}}.geneBodyCoverage_plot.r"
    ],
    sample=SAMPLES
)

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
        expand(
            f"{RSEQC_DIR}/{{sample}}.bam_stat.txt",
            sample=SAMPLES
        ),
        expand(
           f"{RSEQC_DIR}/{{sample}}.read_distribution.txt",
           sample=SAMPLES
        ),
        GENEBODY_OUTPUTS,
        INNER_DISTANCE_OUTPUTS,
        SAMPLE_COMPARISON_OUTPUTS

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

rule bam_stat:
    input:
        bam=f"{RESULTS_DIR}/alignment/{{sample}}.sorted.bam"

    output:
        report=f"{RSEQC_DIR}/{{sample}}.bam_stat.txt"

    log:
        f"{LOGS_DIR}/bam_stat/{{sample}}.log"

    conda:
        "environment/rseqc.yaml"

    shell:
        """
        mkdir -p {RSEQC_DIR} $(dirname {log})

        bam_stat.py \
            -i {input.bam} \
            > {output.report} 2> {log}
        """

rule read_distribution:
    input:
        bam=f"{RESULTS_DIR}/alignment/{{sample}}.sorted.bam",
        bed12=BED12

    output:
        report=f"{RSEQC_DIR}/{{sample}}.read_distribution.txt"

    log:
        f"{LOGS_DIR}/read_distribution/{{sample}}.log"

    conda:
        "environment/rseqc.yaml"

    shell:
        """
        mkdir -p {RSEQC_DIR} $(dirname {log})

        read_distribution.py \
            -i {input.bam} \
            -r {input.bed12} \
            > {output.report} 2> {log}
        """

rule bam_coverage:
    input:
        bam=f"{RESULTS_DIR}/alignment/{{sample}}.sorted.bam"

    output:
        bigwig=f"{BIGWIG_DIR}/{{sample}}.bw"

    log:
        f"{LOGS_DIR}/bam_coverage/{{sample}}.log"

    threads:
        config["threads"]["deeptools"]

    params:
        bin_size=config["quality"]["bigwig_bin_size"],
        normalization=config["quality"]["bigwig_normalization"]

    conda:
        "environment/deeptools.yaml"

    shell:
        """
        mkdir -p {BIGWIG_DIR} $(dirname {log})

        bamCoverage \
            --bam {input.bam} \
            --outFileName {output.bigwig} \
            --outFileFormat bigwig \
            --binSize {params.bin_size} \
            --normalizeUsing {params.normalization} \
            --numberOfProcessors {threads} \
            > {log} 2>&1
        """

rule gene_body_coverage2:
    input:
        bigwig=f"{BIGWIG_DIR}/{{sample}}.bw",
        bed12=BED12_SUBSET

    output:
        report=f"{GENEBODY_DIR}/{{sample}}.geneBodyCoverage.txt",
        pdf=f"{GENEBODY_DIR}/{{sample}}.geneBodyCoverage.pdf",
        rscript=f"{GENEBODY_DIR}/{{sample}}.geneBodyCoverage_plot.r"

    log:
        f"{LOGS_DIR}/gene_body_coverage2/{{sample}}.log"

    params:
        prefix=lambda wildcards: f"{GENEBODY_DIR}/{wildcards.sample}"

    conda:
        "environment/rseqc.yaml"

    shell:
        """
        mkdir -p {GENEBODY_DIR} $(dirname {log})

        geneBody_coverage2.py \
            -i {input.bigwig} \
            -r {input.bed12} \
            -o {params.prefix} \
            -t pdf \
            > {log} 2>&1
        """

rule subset_bed12:
    input:
        bed12=BED12

    output:
        BED12_SUBSET

    params:
        n=config["quality"]["gene_body_subset_size"],
        seed=config["quality"]["gene_body_subset_seed"]

    log:
        f"{LOGS_DIR}/subset_bed12.log"

    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log})

        python scripts/subset_bed12.py \
            --input {input.bed12} \
            --output {output} \
            --n {params.n} \
            --seed {params.seed} \
            > {log} 2>&1
        """

rule inner_distance:
    input:
        bam=f"{RESULTS_DIR}/alignment/{{sample}}.sorted.bam",
        bed12=BED12

    output:
        stats=f"{INNER_DISTANCE_DIR}/{{sample}}.inner_distance.txt",
        freq=f"{INNER_DISTANCE_DIR}/{{sample}}.inner_distance_freq.txt",
        pdf=f"{INNER_DISTANCE_DIR}/{{sample}}.inner_distance_plot.pdf",
        rscript=f"{INNER_DISTANCE_DIR}/{{sample}}.inner_distance_plot.r"

    params:
        prefix=f"{INNER_DISTANCE_DIR}/{{sample}}",
        sample_size=config["quality"]["inner_distance_sample_size"],
        lower=config["quality"]["inner_distance_lower_bound"],
        upper=config["quality"]["inner_distance_upper_bound"],
        step=config["quality"]["inner_distance_step"],
        mapq=config["quality"]["inner_distance_mapq"]

    log:
        f"{LOGS_DIR}/inner_distance/{{sample}}.log"

    conda:
        "environment/rseqc.yaml"

    shell:
        """
        mkdir -p {INNER_DISTANCE_DIR} $(dirname {log})

        inner_distance.py \
            -i {input.bam} \
            -r {input.bed12} \
            -o {params.prefix} \
            -k {params.sample_size} \
            -l {params.lower} \
            -u {params.upper} \
            -s {params.step} \
            -q {params.mapq} \
            > {log} 2>&1
        """

rule multi_bigwig_summary:
    input:
        bigwigs=expand(
            f"{BIGWIG_DIR}/{{sample}}.bw",
            sample=SAMPLES
        )

    output:
        npz=BIGWIG_SUMMARY,
        raw=BIGWIG_RAW

    log:
        f"{LOGS_DIR}/multi_bigwig_summary.log"

    threads:
        config["threads"]["deeptools"]

    params:
        bin_size=config["quality"]["bin_size"]

    conda:
        "environment/deeptools.yaml"

    shell:
        """
        mkdir -p {COMPARISON_DIR} $(dirname {log})

        multiBigwigSummary bins \
            --bwfiles {input.bigwigs} \
            --smartLabels \
            --binSize {params.bin_size} \
            --numberOfProcessors {threads} \
            --outFileName {output.npz} \
            --outRawCounts {output.raw} \
            > {log} 2>&1
        """

rule plot_correlation:
    input:
        matrix=BIGWIG_SUMMARY

    output:
        pdf=CORRELATION_PDF,
        matrix=CORRELATION_MATRIX

    log:
        f"{LOGS_DIR}/plot_correlation.log"

    params:
        method=config["quality"]["correlation_method"]

    conda:
        "environment/deeptools.yaml"

    shell:
        """
        mkdir -p {COMPARISON_DIR} $(dirname {log})

        plotCorrelation \
            --corData {input.matrix} \
            --corMethod {params.method} \
            --whatToPlot heatmap \
            --plotNumbers \
            --skipZeros \
            --plotTitle "Sample correlation" \
            --plotFile {output.pdf} \
            --outFileCorMatrix {output.matrix} \
            > {log} 2>&1
        """

rule plot_pca:
    input:
        matrix=BIGWIG_SUMMARY

    output:
        pdf=PCA_PDF,
        data=PCA_DATA

    log:
        f"{LOGS_DIR}/plot_pca.log"

    params:
        ntop=config["quality"]["pca_ntop"]

    conda:
        "environment/deeptools.yaml"

    shell:
        """
        mkdir -p {COMPARISON_DIR} $(dirname {log})

        plotPCA \
            --corData {input.matrix} \
            --plotFile {output.pdf} \
            --outFileNameData {output.data} \
            --ntop {params.ntop} \
            --addLabels \
            > {log} 2>&1
        """
