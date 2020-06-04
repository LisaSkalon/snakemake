import pandas as pd

SAMPLES_INFO = pd.read_csv("./data_table.tsv", sep = '\t')
SAMPLES_INFO['Sample'] = SAMPLES_INFO.apply(lambda row: row.GSM+'_'+row.Cell+'_'+row.Target, axis=1)
SAMPLES_INFO.set_index(['Sample'], inplace=True)

SAMPLES = SAMPLES_INFO.index.values.tolist()
#PATH = SAMPLES_INFO['File'].tolist()

config = dict(paths = dict(zip(SAMPLES_INFO.index.values,SAMPLES_INFO.File)))

def html_gen(wildcards):
    return ["qc/fastqc/{}.html".format(i) for i in config['paths'].keys()]
def log_gen(wildcards):
    return ["logs/fastqc/{}.log".format(i) for i in config['paths'].keys()]
def zip_gen(wildcards):
    return ["qc/fastqc/{}_fastqc.zip".format(i) for i in config['paths'].keys()]
def multi_gen(wildcards):
    return ['qc/fastqc/{}_fastqc.zip'.format(i) for i in config['paths'].keys()]

#SAMPLES, = glob_wildcards(SAMPLES_INFO.loc[{sample}, 'File'])

rule all:
    input:
        html_gen,
        zip_gen,
        log_gen,
        "qc/multiqc.html"

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
        "qc/fastqc/."
    output:
        "qc/multiqc.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log"
    wrapper:
        "0.59.2/bio/multiqc"