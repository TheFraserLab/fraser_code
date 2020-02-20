# Commonly used code from the Fraser Lab

This repo is intended to be a collection of general-purpose software, scripts, and functions written by members of the Fraser Lab at Stanford University.

The code deposited here is diverse, ranging from bioinformatic tools to generic plotting scripts.

We intend to increase efficiency in our workflows by reducing the amount of time we spend solving common issues we encounter in the dry lab, as it is common that several people have worked on similar coding problems in the past.

We also aim to transparently share our strategies when dealing with common tasks in our lab, e.g. processing beginning-to-end data in the realm of gene expression, allele-specific expression, etc. 

## Folder structure

We organize our code by topics, one per folder with each containing its own readme and likely having specific guidelines for contributing. 

We currently cover the following topics:

- [compbio_snippets](https://github.com/TheFraserLab/fraser_code/tree/master/compbio_code_snippets) — miscellaneous code used to solve common bionformatic issues (e.g. dealing with bam files).
- [config_files](https://github.com/TheFraserLab/fraser_code/tree/master/config_files) — config files for common software, e.g. vim, bash, etc.
- [gtex](https://github.com/TheFraserLab/fraser_code/tree/master/gtex) — code to process and access [GTEx](https://gtexportal.org/home/) data.
- [plots](https://github.com/TheFraserLab/fraser_code/tree/master/plots) — whether it's R or python we have a solution for your visualization task.
- [tcga](https://github.com/TheFraserLab/fraser_code/tree/master/tcga) - code to process and access [TCGA](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga) data.
- [workflows](https://github.com/TheFraserLab/fraser_code/tree/master/workflows) — a collection of our pipelines in snakemake or nextflow.

## Contribute

We ask you to fork the repo to your own account, make changes and then perform a pull request.

**Please read the guidelines below before submitting a pull request**

### Step-by-step process

#### File/script/function

1. Fork this repo to your github account.
2. Clone to your laptop or a server.
3. Create or copy a file in the appropriate folder.
4. Update the lookup list of the corresponding readme.
5. Push to your fork.
6. Submit a pull request to this repo.

#### Folder

1. Fork this repo to your github account.
2. Clone to your laptop or a server.
3. Create a new folder
3. [Optional] Create or copy a file in the new folder.
4. Create a readme and a lookup table for the folder.
5. Push to your fork.
6. Submit a pull request to this repo.

#### Workflow
1. The workflow should be hosted in its own repo.
2. Add the link to that repo in the workflow's readme.

## Guidelines

While we provide general guidelines that all code should follow, each folder may have its specific guidelines, so please make sure to take a review the corresponding readme.

### Content
- One file per issue, and if possible one function per file. *Each file should only contain code to solve a **specifc** task. Multiple function are allowed only if they are all contributing to the same issue. Please avoid submitting a "collection" of. functions*
- Each file and its description has to be included in the lookup list of the corresponding folder's readme.
- Workflows should be in the format of a pipelining software like nextflow or snakemake. We only include links to their corresponding repos (i.e. your pipeline repo).

### Documentation
- All functions/scripts must contain:
 - A brief description of what they do.
 - Explanation of inputs, outputs, and arguments.
 - List of dependencies.
 - An example, if possible.
- Follow these conventions:
 - Python
     - Docstrings ([any style](https://stackoverflow.com/questions/3898572/what-is-the-standard-python-docstring-format?answertab=active#tab-top) is ok).
     - [Type hints](https://docs.python.org/3/library/typing.html) are not required but highly recommended.
 - R
     - [roxygen](https://r-pkgs.org/man.html#roxygen-comments).
 - shell/bash
     - (suggestions?).
- Workflows
  - Include the description along a brief description of important steps included in the pipeline that can be useful to other people (e.g. mapping reads).
  - See the corresponding folder for more details.
