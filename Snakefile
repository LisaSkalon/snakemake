import pandas as pd

configfile: "config.yaml"

#ruleorder: bowtie2 > multiqc_bams

SAMPLES_INFO = pd.read_csv("./data_table.tsv", sep = '\t')
SAMPLES_INFO['Sample'] = SAMPLES_INFO.apply(lambda row: row.GSM+'_'+row.Cell+'_'+row.Target, axis=1)
SAMPLES_INFO.set_index(['Sample'], inplace=True)

SAMPLES = SAMPLES_INFO.index.values.tolist()
GENOME=config["genome"]

def multi_gen(wildcards):
    return ['qc/fastqc/{}_fastqc.zip'.format(i) for i in SAMPLES]
def log_gen(wildcards):
    return ["logs/bowtie2/{}_{}.log".format(i, GENOME) for i in SAMPLES]

rule all:
    input:
        #expand("qc/fastqc/{sample}.html", sample = SAMPLES),
        #expand("qc/fastqc/{sample}_fastqc.zip", sample = SAMPLES),
        #expand("logs/fastqc/{sample}.log", sample = SAMPLES),
        "qc/multiqc/reads.html",
        expand('indexes/{genome}/{genome}.fa.gz', genome=GENOME),
        expand("bowtie2-index/{genome}.1.bt2", genome=GENOME),
        expand("bowtie2-index/{genome}.2.bt2", genome=GENOME),
        expand("bowtie2-index/{genome}.3.bt2", genome=GENOME),
        expand("bowtie2-index/{genome}.4.bt2", genome=GENOME),
        expand("bowtie2-index/{genome}.rev.1.bt2", genome=GENOME),
        expand("bowtie2-index/{genome}.rev.2.bt2", genome=GENOME),
        expand("bams/{sample}_{genome}.bam", sample = SAMPLES, genome = GENOME),
        "qc/multiqc/bams.html",
        expand("logs/bowtie2/{sample}_{genome}.log", sample = SAMPLES, genome=GENOME),
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
    threads: 4
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
    priority: 100
    shell:
        #!!!! if you want to use full genome (hg19, ...), uncomment line below !!!!
       #"""wget http://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.genome}/bigZips/{wildcards.genome}.fa.gz \
       #   -O {output}"""
       #!!!! if you want to use full genome (hg19, ...), comment line below !!!!
       """wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/{wildcards.genome}.fa.gz  \
       -O {output}"""

#rule index_bowtie:
#    output:
#        directory('bowtie2-index')
#    params:
#        files_list='indexes/{genome}/{genome}.fa.gz'.format(genome=config['genome']),
#        target='bowtie2-index/{genome}'.format(genome=config['genome'])
#    shell:
#        """mkdir -p {output} $$ bowtie2-build {params.files_list} {params.target}"""

rule index_bowtie:
    input:
        rules.wget.output
    output:
        output1="bowtie2-index/{genome}.1.bt2",
        output2="bowtie2-index/{genome}.2.bt2",
        output3="bowtie2-index/{genome}.3.bt2",
        output4="bowtie2-index/{genome}.4.bt2",
        outputrev1="bowtie2-index/{genome}.rev.1.bt2",
        outputrev2="bowtie2-index/{genome}.rev.2.bt2"
    priority: 50
    threads: 4
    params:
        basename="bowtie2-index/{genome}"
    shell: "bowtie2-build {input} {params.basename}"

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
   priority: 25
   wrapper:
       "0.60.0/bio/bowtie2/align"

# rule bowtie2:
#     input:
#         lambda wildcards:
#            expand("reads/{files}", files = SAMPLES_INFO.loc[wildcards.sample, 'File'], sample=SAMPLES)
#     output:
#         "bams/{sample}_{genome}.bam"
#     params:
#        index=rules.index_bowtie.params.basename
#     shell:
#         "bowtie2 -x {params.index} -p 4 -U {input} -S {output}"


rule multiqc_bams:
   input:
       log_gen
   output:
         dir="qc/multiqc",
         name="bams.html"
  # !!! environment: environment.yaml : you don`t need specific conda env for this rule,
  # !!! just the one already mentioned
   priority: 10
   shell:
       """multiqc --force -o {output.dir} -n {output.name} {input}"""

rule samtools_sort:
    input:
         "bams/{sample}_{genome}.bam"
    output:
         "bams/sorted/{sample}_{genome}.sorted.bam"
    params:
        ""
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@.
    priority: 8
    wrapper:
        "0.60.0/bio/samtools/sort"

rule samtools_index:
    input:
        "bams/sorted/{sample}_{genome}.sorted.bam"
    output:
        "bams/sorted/{sample}_{genome}.sorted.bam.bai"
    params:
        "" # optional params string
    priority: 7
    wrapper:
        "0.60.0/bio/samtools/index"

rule bigWig:
    input:
        sorted = "bams/sorted/{sample}_{genome}.sorted.bam",
        index ="bams/sorted/{sample}_{genome}.sorted.bam.bai"
    output:
        dir="bigwig/{sample}_{genome}.bw"
    threads: 4
    shell:
        """
        bamCoverage \
        --bam {input.sorted} \
        --outFileName {output.dir} \
        --binSize 1 \
        --outFileFormat bigwig \
        """

# rule tar_gz:
#     input:
#         bigWig = rules.bigWig.output,
#         multiqc_reads = rules.multiqc_bams.output.dir
#         multiqc_bams = rules.multiqc_bams.output
#     output:
#         "chip-seq_{genome}.tar.gz"
#     shell:
#         """tar -czvf {output} {input.multiqc_reads} {input.bigWig} """
