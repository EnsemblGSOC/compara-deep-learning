# Gene Pair Orthology Prediction Using Deep Learning

## Project Goal

The aim of this project is to take pairs of genes from different species and predict the nature of their evolutionary relationship. Primarily, this concerns whether the genes are orthologs or paralogs (see: [here](https://useast.ensembl.org/info/genome/compara/homology_types.html)). The compara project at the EBI currently uses a method matches sequence similarity and the species tree to categorise gene pairs. Whilst it has proven successful in the past, the current methods do not scale as the number of comparison pairs increases. Thus, as the genomes of more species, the ability to assess the evlutionary relationship between genes is inhibited by the computational resources required for this task. 

This project is a proof of principle that <em>deep learning</em> is capable of predicting orthology relationships with a high degree of accuracy between even distant species, with a view to demonstrating it's viability for predicting orthology at scale in a production setting

## The pipeline

The pipeline ![alt text](https://github.com/AidanMar/compara-deep-learning/blob/master/pipeline3/small_dag.png)

