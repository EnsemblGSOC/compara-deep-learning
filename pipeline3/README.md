# GSoC (2021) - Gene Pair Orthology Prediction Using Deep Learning

## Project Goal

The aim of this project is to take pairs of genes from different species and predict the nature of their evolutionary relationship. Primarily, this concerns whether the genes are orthologs or paralogs (see: [here](https://useast.ensembl.org/info/genome/compara/homology_types.html)). The compara project at the EBI currently uses a method matches sequence similarity and the species tree to categorise gene pairs. Whilst it has proven successful in the past, the current methods do not scale as the number of comparison pairs increases. Thus, as the genomes of more species, the ability to assess the evlutionary relationship between genes is inhibited by the computational resources required for this task. 

This project is a proof of principle that <em>deep learning</em> is capable of predicting orthology relationships with a high degree of accuracy between even distant species, with a view to demonstrating it's viability for predicting orthology at scale in a production setting. 

## The network's input

In the project, a neural network capable of attaining > 99% validation set accuracy was attained using a <em>convolutional neural network</em> which used a variety of features in the form of matrices as well as the distance between species to form the basis of it's predictions.

Every input matrix represents the genes in the neighborhood of the two genes being compared. Let's call the comparison genes, whose orthology we wish to assess, A and B. In the current implementation the 3 genes upstream and the 3 genes downstream of both genes A & B on their respective chromosomes are obtained, giving two vectors of 7 genes. For instance, the vector representing gene A and its neighbours can be written as v1 = [d3,d2,d1,A,u1,u2,u3]. Here the d elements represent the downstream genes, A represents gene A itself, and the u elements represent the upstream genes. v2 is likewise formed for gene B. All pairwise combinations of genes between these two vectors are then made, and then stored in a 7x7 matrix, called M. As such, element i,j of M represent some comparison between gene i from V1 and gene j from V2. Several different comparison mertics are used, resulting in multiple matrices, each housing different metrics of comparison between the genes in the neighbourhood of genes A and B. 5 such types matrices are formed which are meant to capture evolutionary information relevant to the development of classifying the orthology relationship. These feature matrices are briefly as follows:

1) The global alignment matrix
2) The local alignment score matrix
3) The local alignment coverage matrix
4) The Pfam Jaccard matrix
5) The Pfam coverage matrix

Further details of each of these matrices can be found in the steps detailing the pipelines operation.

## The network

For each type of these feature matrices, the network has a separate CNN sub-module that learns a feature map for that matrix. (Each feature get's its own CNN because the features are so diverse there's no reason to for them to apply the same filters to each of these layers. All features are later joined in MLP layers to merge all of them together). Each CNN sub-module is visualised immediately below. The output of the sub-module is a fully connected MLP layer derived from the earlier CNN matrix values. 

![alt text](https://github.com/AidanMar/compara-deep-learning/blob/master/pipeline3/CNN_model.png)

The output for each sub-module is then concatanated to and fed into an MLP layer is seen below, which in turn is concatenated with the value of the species distance between the species containing the two comparison genes. All of these are then either fed into a two or three way softmax layer which is used to assign probabilities for the gene pairs orthology category. The two way classification could simply be ortholog vs paralog, whilst a three way classification task would be ortholog vs paralog vs non-homologous. The nodes in the first layer of the network are coloured by their feature type, and by extension, the CNN sub-module that they originated from (with only 5 input nodes representing the 50 for cleaner visualisation).

![alt text](https://github.com/AidanMar/compara-deep-learning/blob/master/pipeline3/MLP_layers.png)


## Snakemake  implementation

The project is implemented in [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) and makes use of a conda environment stored in a [yml](https://github.com/AidanMar/compara-deep-learning/blob/master/pipeline3/compara.yml) in this repo to ensure as much reproducibility as possible. Moreover, once this conda environment is installed, all packages required for running the pipeline from end to end should already be installed.

Snakemake allows us to specify a sequence of <em>rules</em> that can be used to transform input files into output files. By stringing many such rules together, snakemake allows us to elegantly formulate an entire data processing pipeline. Whilst the sequence of rules and transformations can be challenging to understand by reading a pipeline's script, the Directed Acyclic Graph of jobs that are used to go from raw input to output can be visualised to aid in in understanding the pipelines steps. An example of this can be seen here:

![alt text](https://github.com/AidanMar/compara-deep-learning/blob/master/pipeline3/small_dag.png)

At first, this can be a little intimidating to view, however after a little thought it can be seen how much this visualisation aids in understanding the process by informing us which rules in the pipeline are dependent on which others before they run. An arrow from one rule to the next indicates the sequence in which a rule must be run. Moreover, once a rule has been established, it can further be generalised across multiple examples. For instance, above, the DAG is demonstrated for <em>neolamprologus brichardi</em>. But we can easily run the pipeline for multiple species, for example, the dag below generalises this to 3 species. We can in principle generalise the pipeline to as many species as we would like, provided that we have the necessary raw data and compute power to execute the pipeline. For the current demonstration, the pipleine was run with 56 species(the DAG can be seen [here](https://github.com/AidanMar/compara-deep-learning/blob/master/pipeline3/dag.png), though the image is very large).

![alt text](https://github.com/AidanMar/compara-deep-learning/blob/master/pipeline3/medium_dag.png)

## The pipeline

### Pre-requisites
If you pull the repo, onto the EBI's LSF cluster system, called Codon, it should be able to run after some initial set-up. 

To get going, make sure you have access to anaconda3 or miniconda so that you can work with conda environments. Firstly, clone this git repo to a suitable place on your cluster system. Once this has been complete, install the conda environment specified by the compara.yml in this repo. To create the environment, simply run:
```
 conda env create -f compara.yml
```
There are a lot of packages requires to run this pipeline from end to end, so do be patient with the installation process. For more info on dealing with conda envs and yml files see [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file).

### Setting up Diamond

If you have installed the conda environment, you have installed [Diamond](https://github.com/bbuchfink/diamond). This is the tool used for align genes to the PFAM-A database and is used in the pipeline to generate a database of alignments for each species. Before this can be done however, Diamond require some pre-setup to get going. Go to the ebi [ftp](http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/) and download Pfam-A.fasta.gz. Then run 
```
diamond makedb --in /path1/Pfam-A.fasta.gz -d /path2/Pfam-A
```
where the inputs and outputs are modifiable to whereever is most appropriate on your system.

### Config File

Now that the package environment is set-up we need to specify the config file for the environment. This is the config.json file. By modifying this file, most of the flow of information through the pipeline should be readily controlled.

- species_file: contains the list of species that we wish to sample homologous pairs from
- reports_out_base_path: The path to the directory which will then contain all the subdirectories to which a stdout will write for each of the jobs that snakemake submits to the cluster. For instance, the rule split_diamond_alignment will will submit jobs to the cluster and then write its stdout to config["reports_out_base_path"] + "/Diamond_Split/"
- num_neighbours: Controls the number of upstream and downstream genes used in generated the feature matrices. This will results (2n + 1) x (2n + 1) feature matrices.
- species_neighbour_outdir: Path to directory of where to store neighbour genes as a dictionary for each species
- Longest_pep_fasta_path: Path to directory of where to store the longest peptide sequences of every gene
- Diamond_alignments_path: Path to directory of where to store the output of the Diamond Alignments
- Diamond_Pfam_db_path: Path to the Pfam-A database generated in the Diamond configuration step above. Ie, wherever you put /path2/Pfam-A earlier.
-samples_neighbour_outdir: Path to directory containing all the gene which are randomly sampled from the homology databases
- samples_per_species: The number of samples you wish to get from each species homology database
- num_homolog_species: The number of homology species to attempt to get from each database
- sp_names: Path to the file containing names of the species in the species distance matrix. This is included in this repo
-dist_matrix: Path to the file containing the species distance matrix. This is included in this repo
- out_dir: Path to directory where you would like most of the output information to go to.
- Final_data: Path to the directory where the final compiled and concatednated will be save. These files are ~1-10Gbs in size.
- Model_outputs: Path to the directory where the saved models will be saved as well as their evaluation statistics.

### Config directory

Within this repo there is a directory called config. This has a set of files, each one containing the file paths to files needed to run this pipeline. This was run natively on the EBI server, and therefore all files are readily available on the system. If you're running this on another system, then you'll need to download these files to wherever you run this pipeline.

The three main file types are the <em>peptide</em> files which contain the protein sequences, the <em>CDS</em> files which hold the coordinates of the protein coding sequences of each genes, and the <em>homology database</em> files which contain the list of labelled homologous pairs for each species that the EBI has data from. In order to run this pipeline for a species, you need to have all three such files for your species of interest. If not, it cannot be run. By investigating the species which have all three files present, you can compile a list of species you wish to run the pipeline for. I did this, and placed a completed list into config/valid_species.txt. My selection procedure was relatively random and I just picked ~50 species that had all the relevant data. You can generate your own species set if you're interested in other species.

### Running the Pipeline

Everything in the pipeline is controlled by the snakefile. Each rule of the snakefile can be seen earlier in the repo. Here is described what each rule does, in their order in the DAG, and what script in the scripts/ folder is required for each rule. Whilst snakemake implements each of these rules in turn, it simply strings together a sequence of bash and python scripts which are modularised such that one can in principle take the script and run each step separately without snakemake if desired. Some modification to the input output format for each script would be required to enable this to happen.

#### rule: select_longest

Inspects each species fasta file, selects the longest version of each gene and writes it to a knew fasta file.

script: GeneReferencePetpideSelection.py

#### rule: Diamond_alignment

Takes each of the fasta files output from the rule<em>select longest</em> and aligns each genes for that species to the Pfam-A database to identify all of the Pfam domains within each gene. The output is a tsv file containing each gene, and the pfam domain that it aligns to.

script: uses diamond which is part of the compara.yml conda environment

rule: split_diamond_alignment

Each of the tsv files output from a diamond alignment can be very large, too large to hold in memory. Therefore, each tsv file is split into many separate tsv files, each containing the pfam alignments for each gene, as these can easily be loaded into memory one at a time.

#### rule: species_neighbour

Goes through every gene in a species genome using the species' gtf file, and finds its upstream and downstream neighbouring genes. These are then written to a dictionary. This is done serparately for each species.

script: geneNeighbourFormater.py

#### rule: select_samples

Reads the homology tsv file of the species, and selects homologous gene pairs between that species the homology file belongs to and multiple other "query species". There are many possible ways to select candidate gene pairs from this tsv file. The current implementation reads in the valid_species.txt file, and only allows comparison between species in that list. The potential query species are then ranked in order of their species distance, as stated in the species distance matrix. N query species are then chosen at evenly spaced ranks, where n is controlled in the config.json by "num_homolog_species". The rule aims to sample roughly equal numbers of orthologs and paralogs for comparison. Non-homologous genes are also randomly sampled between this set of species to give a list of gene pairs and their orthology type labels of ["ortholog","paralog","not"]. This is not the only way to generate a candidate genes pairs dataset, but the random sampling should prevent forms of bias slipping into the dataset. The script mostly exploits pandas to do the heavy lifting. The output pair names, along with some other features about the pairs, are saved as csv files at the location specified by "samples_neighbour_outdir". Additionally, each sampled gene's neighbourhood is identified by querying the genes name against the dictionaries produced by the rule <em>species_neighbour</em>. These are saved in {SPECIES}_sample_neighbour.json as dictionaries. Similarly, this is done for the each gene's homolog {SPECIES}_homology_sample_neighbour.json. The "negative" examples are those that are not homologous pairs, but just randomly sampled, unrelated pairs of genes. These are saved at {SPECIES}_negative_sample_neighbour.json, and {SPECIES}_negative_homology_sample_neighbour.json respectively. 

script: sample_gene.py

#### rule: synteny_matrix & negative_synteny_matrix

This rule produced the first set of feature matrices that will in turn be fed into the network. The index of each element in the matrix was described earlier in this README under the heading <em>The network's input</em>. We use the same notation here. These matrices saved with the following names as numpy's npy files:

{SPECIES}_global1.npy
{SPECIES}_global2.npy
{SPECIES}_local1.npy
{SPECIES}_local2.npy

Gene comparisons of global1 are made by taking the normalised levenstein matrix distance between each gene i and gene j. Matrix global2 is identical, except the sequences of gene j is now reversed. 

local1 in fact contains 3 matrices, stacked one on top of another like the RGB channels of an image. The first of these contains a normalised Smith-Waterman alignment score between the compared genes. The next matrix contains the normalised coverage of compared gene i against compared gene j. The final matrix contains the normalised coverage of gene j against gene i as the coverage opration is not symmetric. local2 is identical to local1 except that gene j is has its sequence reversed.

The "negative" versions of these files are identical to those described above, except it is ran for those pairs of genes that were identified as being non-homologous/unrelated in the <em>select_samples</em> step.

script: synteny_matrix.py & negative_synteny_matrix.py

#### rule: pfam_matrix & negative_pfam_matrix

This set of matrices utilises the Diamond alignments to the Pfam-A, and therefore requires both sample genes and the sequence of steps relating to diamond to have been run. This is reflected in the earlier DAG. 

The output is saved as {SPECIES}_negative_Pfam.npy. This contains 3 matrices stacked like an RGB channel. These Pfam matrices are aimed to contain information regarding the kind of functional structural domains that are shared between two genes, which should enable orthologs and paralogs to be differentiated from one another. Earlier in the pipeline, we found the list of Pfam domains that each gene aligned, and the positions of those alignments, using Diamond. When comparing gene i and gene j at this step, take the set of Pfam-domains each gene aligns to, along with the coordinates of those alignments. We then compute the [Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index) between these sets where in order for a Pfam gene to be in both sets, both genes must have aligned to it, and at an overlapping position. A Jaccard index of 1 indicates the genes share lot of functional domains, whilst a Jaccard index of 0 indicates no such shared domains are present, and this should be reflective of the different natures of the evolitionary relationships between different gene pairs. This information forms the first matrix in the pfam.npy.

Whilst likely useful, the Jaccard index as a scalar value is being asked to represent a lot of degrees of variation in the range of functional domains shared between evolutionary relationships. The next two matrices in the pfam.npy aim to rectify this by taking the coverage overlap between the two pairs of sequences, which aims to capture variable information that might reflect, for instance, disparities in the lengths of the two genes. Again the coverage of A on B is not the same as B on A, so two matrices are requried to represent this in both directions.

Once more, this is done separately for the "negative"/non-homologous pairs, by using a different script.

script: Pfam_matrix_parallell_split.py & negative_Pfam_matrix_parallel_split.py



## So what's in each of the matrices?
