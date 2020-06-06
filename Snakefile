import pandas as pd

configfile: "config.yaml"

SAMPLES_INFO = pd.read_csv("./data_table.tsv", sep = '\t')
SAMPLES_INFO['Sample'] = SAMPLES_INFO.apply(lambda row: row.GSM+'_'+row.Cell+'_'+row.Target, axis=1)
SAMPLES_INFO.set_index(['Sample'], inplace=True)

SAMPLES = SAMPLES_INFO.index.values.tolist()
GENOME=config["genome"]

def multi_gen(wildcards):
    return ['qc/fastqc/{}_fastqc.zip'.format(i) for i in SAMPLES]

rule all:
    input:
        expand("qc/fastqc/{sample}.html", sample = SAMPLES),
        expand("qc/fastqc/{sample}_fastqc.zip", sample = SAMPLES),
        expand("logs/fastqc/{sample}.log", sample = SAMPLES),
        "qc/multiqc/reads.html",
        expand('indexes/{genome}/{genome}.fa.gz', genome=GENOME),
        expand("bams/{sample}_{genome}.bam", sample = SAMPLES, genome = GENOME),
        "qc/multiqc/bams.html",
        expand("bams/sorted/{sample}_{genome}.sorted.bam", sample = SAMPLES, genome = GENOME),
        expand("bams/sorted/{sample}_{genome}.sorted.bam.bai", sample = SAMPLES, genome = GENOME),
        expand("bigwig/{sample}_{genome}.bw", sample = SAMPLES, genome = GENOME),
        expand('chip-seq_{genome}.tar.gz', genome = GENOME)

rule fastqc:
    input:
        lambda wildcards:
        expand("reads/{files}", files = SAMPLES_INFO.loc[wildcards.sample, 'File'], sample=SAMPLES)
    output:
        html='qc/fastqc/{sample}.html',
        zip='qc/fastqc/{sample}_fastqc.zip',
    params: ""
    log:
        'logs/fastqc/{sample}.log'
    wrapper:
        "0.59.2/bio/fastqc"


rule multiqc_reads:
    input:
        multi_gen
    output:
        "qc/multiqc/reads.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log"
    wrapper:
        "0.59.2/bio/multiqc"

rule wget:
    output:
        'indexes/{genome}/{genome}.fa.gz'
    shell:
        #!!!! if you want to use full genome (hg19, ...), uncomment line below !!!!
       #"""wget http://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.genome}/bigZips/{wildcards.genome}.fa.gz \
       #   -O {output}"""
       #!!!! if you want to use full genome (hg19, ...), comment line below !!!!
       """wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/{wildcards.genome}.fa.gz  \
       -O {output}"""

rule index_bowtie:
    output:
        directory('bowtie2-index')
    params:
        files_list='indexes/{genome}/{genome}.fa.gz'.format(genome=config['genome']),
        target='bowtie2-index/{genome}'.format(genome=config['genome'])
    shell:
        """mkdir -p {output} $$ bowtie2-build {params.files_list} {params.target}"""

rule bowtie2:
    input:
         sample=
         lambda wildcards:
            expand("reads/{files}", files = SAMPLES_INFO.loc[wildcards.sample, 'File'], sample=SAMPLES)
    output:
        "bams/{sample}_{genome}.bam"
    log:
        "logs/bowtie2/{sample}_{genome}.log"
    params:
        index="bowtie2-index/{genome}",  # prefix of reference genome index (built with bowtie2-build)
        extra=""  # optional parameters
    threads: 4  # Use at least two threads
    wrapper:
        "0.60.0/bio/bowtie2/align"

rule multiqc_bams:
    input:
        "logs/bowtie2/"
    output:
        "qc/multiqc/bams.html"
   #! set environment conda:
    shell:
        """multiqc --force -n {output} {input}"""

rule samtools_sort:
    input:
         "bams/{sample}_{genome}.bam"
    output:
         "bams/sorted/{sample}_{genome}.sorted.bam"
    params:
        ""
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@.
    wrapper:
        "0.60.0/bio/samtools/sort"

rule samtools_index:
    input:
        "bams/sorted/{sample}_{genome}.sorted.bam"
    output:
        "bams/sorted/{sample}_{genome}.sorted.bam.bai"
    params:
        "" # optional params string
    wrapper:
        "0.60.0/bio/samtools/index"

rule bigWig:
    input:
        sorted = "bams/sorted/{sample}_{genome}.sorted.bam",
        index ="bams/sorted/{sample}_{genome}.sorted.bam.bai"
    output:
        "bigwig/{sample}_{genome}.bw"
    shell:
        """
        bamCoverage \
        --bam {input.sorted} \
        --outFileName {output} \
        --binSize 1 \
        --outFileFormat bigwig \
        """

rule tar_gz:
    input:
        bigWig = "bigwig/",
        multiqc = "qc/multiqc"
    output:
        "chip-seq_{genome}.tar.gz"
    shell:
        """tar -czvf {output} {input.bigWig} {input.multiqc} """
