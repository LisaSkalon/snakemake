# ChIP-Seq processing mini-pipeline (Snakemake)

### Configuration

The only tool required to launch the pipeline is `conda`.  

    If conda is not installed, follow the instructions at Conda website.
    Navigate to repository directory.

Create a Conda environment for snakemake:  

`$ conda env create --file ./environment.yaml --name snakemake`

Activate the newly created environment:

`$ source activate snakemake`

### Launch

To launch the pipeline, type the line below from your working directory, there you store sample info file

`snakemake --cores 1 --use-conda`

Genome=chr15 (from hg19) by default. If you want to change reference genome, type `--config genome=<genome>` , <genome> = desired genome. Uncomment line 72 in Snakemake file (see detailed instructions in the file)
  
### Output

Archive with output (multiqc reports & bigwig files): https://drive.google.com/drive/folders/1Lr3sQfBaDoOZTHgZwy1atA_9qGCLldxo?usp=sharing


### Author
Elizaveta Skalon
