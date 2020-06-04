import pandas as pd

SAMPLES_INFO = pd.read_csv("./data_table.tsv", sep = '\t')
SAMPLES_INFO['Sample'] = SAMPLES_INFO.apply(lambda row: row.GSM+'_'+row.Cell+'_'+row.Target, axis=1)

config = dict(paths = dict(zip(SAMPLES_INFO.Sample,SAMPLES_INFO.File)))

def html_gen(wildcards):
    return ["qc/fastqc/{}.html".format(i) for i in config['paths'].keys()]
def log_gen(wildcards):
    return ["logs/fastqc/{}.log".format(i) for i in config['paths'].keys()]
def zip_gen(wildcards):
    return ["qc/fastqc/{}_fastqc.zip".format(i) for i in config['paths'].keys()]


rule all:
    input:
        html_gen,
        zip_gen,
        log_gen

rule fastqc:
    input:
        expand('reads/{path}', path = config['paths'].values())
    output:
        html='qc/fastqc/{sample}.html',
        zip='qc/fastqc/{sample}_fastqc.zip',
    params: ""
    log:
        'logs/fastqc/{sample}.log'
    wrapper:
        "0.59.2/bio/fastqc"
