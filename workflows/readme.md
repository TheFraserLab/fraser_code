# Fraser lab workflows

This a collection of workflows used in the Fraser lab for complex bioinformatic processes.

Specific steps of a workflow can be useful in other settings, e.g. mapping reads to a genome, calculating gene differential expression, GO enrichments, etc.

## List of workflows

- **[Create plink bed files from VCF files](https://github.com/pablo-gar/create_plink_bed_1000Genomes)**
	- **Description**: From VCF files created from different people, this creates bed files that are used by plink as inputs for most of its commands. 
	- **Pipeline**: snakemake
	- **Author**: Pablo Garcia
	- **Commonly used features**: NA

## Contribute

If you wish to contribute, please modify this readme and add a link to the repo of your workflow

Please follow these guidelines:

- Only links are provided, we are not hosting the workflows here.
- Please keep links order alphabetically.
- We strongly recommend using a pipelining software like snakemake or nextflow.
- Provide a brief description, the author, the pipelining software, and commonly used features that can be useful outside the pipeline (e.g. mapping reads step)
