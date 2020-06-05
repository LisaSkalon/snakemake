import pandas as pd

SAMPLES_INFO = pd.read_csv("./data_table.tsv", sep = '\t')
SAMPLES_INFO['Sample'] = SAMPLES_INFO.apply(lambda row: row.GSM+'_'+row.Cell+'_'+row.Target, axis=1)
SAMPLES_INFO.set_index(['Sample'], inplace=True)

SAMPLES = SAMPLES_INFO.index.values.tolist()

config = dict(
    paths = dict(zip(SAMPLES_INFO.index.values,SAMPLES_INFO.File)),
    genome = ['chr15'] )

GENOME=config["genome"]

def html_gen(wildcards):
    return ["qc/fastqc/{}.html".format(i) for i in config['paths'].keys()]
def log_gen(wildcards):
    return ["logs/fastqc/{}.log".format(i) for i in config['paths'].keys()]
def zip_gen(wildcards):
    return ["qc/fastqc/{}_fastqc.zip".format(i) for i in config['paths'].keys()]
def multi_gen(wildcards):
    return ['qc/fastqc/{}_fastqc.zip'.format(i) for i in config['paths'].keys()]


ruleorder: bowtie2 > multiqc_bams

rule all:
    input:
        html_gen,
        zip_gen,
        log_gen,
        "qc/multiqc.html",
        expand('indexes/{genome}/{genome}.fa.gz', genome=GENOME),
        expand("bams/{sample}_{genome}.bam", sample = SAMPLES, genome = GENOME),
        "qc/multiqc/bams.html"

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


rule multiqc:
    input:
        multi_gen
    output:
        "qc/multiqc.html"
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
       # """wget http://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.genome}/bigZips/{wildcards.genome}.fa.gz -O {output}""""
       """wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/{wildcards.genome}.fa.gz -O {output}"""

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
        """multiqc -n {output} {input}"""


