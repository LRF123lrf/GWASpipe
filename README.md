# GWASpipe
An automated GWAS analysis toolkit.


# 1.Installation
====================

   
1.1 Install Miniconda
As this workflow is based on the workflow management system snakemake and conda. We strongly recommend installing miniconda3 with python3. For specific installation methods and usage methods of GWASpipe, please refer to the GWASpipe_Documentation.pdf.

1.2 Install GWASpipe
Clone the repository:

```bash
git clone https://github.com/LRF123lrf/GWASpipe.git
```




Create the environment:

```bash
conda env create -n GWASpipe -f envs/envs.yaml
```

Activate the environment:

```bash
conda activate GWASpipe
```

# 2.Prepare input files
====================

Several input files are required in order to run the workflow, a genome file (.fa), an  phenotype file (.csv) and compressed sequencing files (.fastq.gz)

| File type    | Description                                                         |
|--------------|---------------------------------------------------------------------|
| genome.fa    | user-provided genome file containing the genome sequence            |
| .fastq.gz    | user-provided compressed sequencing files                            |
| phenotype.csv    | user-provided compressed phenotype data                            |
| config.yaml  | configuration file to customize the workflow                        |


2.1 Annotation.gff and genome.fa
We recommend retrieving both the genome and the annotation files for your organism from National Center for Biotechnology Information (NCBI) Note: if you use custom annotation files, ensure that you adhere to the gtf standard

2.2 Input .fastq files
These are the input files provided by the user. Both single end and paired end data is supported. Note: Please ensure that you compress your files in .gz format and .fastq.gz

2.3 Set up configuration
Modify the metafile describing your data configs/phenotype.csv .

| Sample ID | Phenotype             |
|-----------|-----------------------|
| S0001     | 0.127029604293491     |
| S0002     | NA                    |
| S0003     | NA                    |
| S0004     | -0.289023784002341    |
| S0005     | 1.02183456721831      |
| S0006     | NA                    |

Customize the workflow based on your need in configs/config.yaml .It contains the following variables:

- **PROJECT**: Project name
- **READSPATH**: The path to fastq files
- **SAMPLES**: configs/phenotype.csv
- **END**: sequencing paired-end or single-end
- **OUTPUTPATH**: The path for final outputs
- **GENOME**: The path of genome files
- **ANNOTATION**: The path of annotation files


# 3. Run RNApipe
====================

   
```bash
sh snakemake.sh
```
