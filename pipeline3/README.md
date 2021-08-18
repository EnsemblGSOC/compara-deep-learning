# Gene Pair Orthology Prediction Using Deep Learning

## Project Goal

The aim of this project is to take pairs of genes from different species and predict the nature of their evolutionary relationship. Primarily, this concerns whether the genes are orthologs or paralogs (see: [here](https://useast.ensembl.org/info/genome/compara/homology_types.html)). The compara project at the EBI currently uses a method matches sequence similarity and the species tree to categorise gene pairs. Whilst it has proven successful in the past, the current methods do not scale as the number of comparison pairs increases. Thus, as the genomes of more species, the ability to assess the evlutionary relationship between genes is inhibited by the computational resources required for this task. 

This project is a proof of principle that <em>deep learning</em> is capable of predicting orthology relationships with a high degree of accuracy between even distant species, with a view to demonstrating it's viability for predicting orthology at scale in a production setting

## The pipeline

The project is implemented in [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) and makes use of a conda environment stored in a [yml](https://github.com/AidanMar/compara-deep-learning/blob/master/pipeline3/compara.yml) in this repo to ensure as much reproducibility as possible. Moreover, once this conda environment is installed, all packages required for running the pipeline from end to end should already be installed.

Snakemake allows us to specify a sequence of <em>rules</em> that can be used to transform input files into output files. By stringing many such rules together, snakemake allows us to elegantly formulate an entire data processing pipeline. Whilst the sequence of rules and transformations can be challenging to understand by reading a pipeline's script, the Directed Acyclic Graph of jobs that are used to go from raw input to output can be visualised to aid in in understanding the pipelines steps. An example of this can be seen here:

![alt text](https://github.com/AidanMar/compara-deep-learning/blob/master/pipeline3/small_dag.png)

At first, this can be a little intimidating to view, however after a little thought it can be seen how much this visualisation aids in understanding the process by informing us which rules in the pipeline are dependent on which others before they run. An arrow from one rule to the next indicates the sequence in which a rule must be run. Moreover, once a rule has been established, it can further be generalised across multiple examples. For instance, above, the DAG is demonstrated for <em>neolamprologus brichardi</em>. But we can easily run the pipeline for multiple species, for example, the dag below generalises this to 3 species. We can in principle generalise the pipeline to as many species as we would like, provided that we have the necessary raw data and compute power to execute the pipeline.

![alt text](https://github.com/AidanMar/compara-deep-learning/blob/master/pipeline3/medium_dag.png)
