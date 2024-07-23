# Pathogen Agnostic Pipeline

This archive contains a viral processing pipeline to generate calls and consensus sequence on any viral pathogen. Designed for viral pathogens this has also been tested on bacteria up to 4MB. 

The pipeline is fully containerised and *should* require only snakemake and apptainer to be installed to run. 

We recomend that this archive is used as a skeleton for further development of pathogen-specific pipelines. To begin this, first generate your new github, then either: 

a) clone this archive, copy all files to your new github, and commit
or
b) add this archive as a second remote and perform a one-way pull, as follows:
  git remote add upstream git@github.com:sethnr/pathag_pipeline.git
  git pull upstream main --allow-unrelated-histories
  git pull upstream main
