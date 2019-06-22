# New Readme

The aim of this project is to use **Deep Neural-Nets** to predict homology type between the give pair of genes. 
This file provides the instructions to replicate the results of this project from scratch.

## Requirements:
1. A machine with atleast **8GB of RAM** (although **16-32GB** is recommended. It depends on the no. of homology databases that you are willing to use in the preparation of the dataset), a graphic card for training the deep neural nets. A single GPU machine would suffice. The model can be trained on CPU as well but will be a lot faster if trained on a GPU.<br/>
2. A stable Internet Connection.<br/>
3. A native/virtual python environment. Install the dependencies from `requirements.txt` using :<br/>
 `pip install -r requirements.txt`<br/>

## Step 1: Data Preparation:
The model uses a synteny matrix and some other factors derived from the species tree to make predictions. <br/>
**Download The Required Files:**
So as to prepare data we will need the following files:<br/>
1.All the `gtf` files to get the start and end locations of the genes and find their neighboring genes. This is used to create the synteny matrix which helps to see the conserved synteny among the genes.<br/>
2. All the `cds` files in `FAST-A` format. They are required but are not mandatory, if the files are not provided the sequences are directly accessed from the `REST API` but the process can be slow :( . It's better to have all the `cds` files.<br/>
3. Homology Databases of Your Choice. All the databases have the same name so its better to change the name of the files with their respective speicies names. One format that works best is `species_name.tsv.gz`.<br/>

To download the `gtf` and `cds` files `ftpg.py` can be used. This scripts writes the links of all the required `gtf` and `cds` files to `gtf-link.txt` and `seq_link.txt`. You can use your own script to download the files or just enter `y` when prompted for permission to download the files. It will automatically download all the files and store them in designated folder. You can manually download each files by pasting link from the files in the browser.<br/>

The designated folders to store the files are as follows:<br/>
`data` all the `gtf` files.<br/>
`data_homology` all the homology databases that you want to use.<br/>
`geneseq`To store the `cds` sequence files in `fasta` format.<br/>
You can assign different folder names if you want but it's better to stick to them. <br/>
**Create Genome Maps:**<br/>
The purpose is to create maps of all the genes present in the `gtf` files with respect to their chromosomes, a map of all the genes belonging to the same chromosome in the given species, a map of all the genes in the given species, a map of all the species whose data has been successfully read.<br/>
To create genome maps run this command:<br/>
`python create_genome_maps.py -d path -r` where:<br/>
`path`: link to the directory where all the `gtf` files exist. If you have used `ftpg.py` then the path is `data`. <br/>
Note: Genome Maps can be downloaded from this [link](https://drive.google.com/open?id=1GjV6dT-Hpf2LWQ-vSpekqqQ7RF_tH8So).<br/>
**Create/Update Neighbor Genes File:**<br/>
This project uses the measure of conserved synteny to predict the homology type. Therefore, to predict the homology type we need to find the neighboring genes of the given homologous pair of genes. <br/>
This file finds the neighboring genes of all the rows in the given databases and writes it to a file called `processed/neighbor_genes.json` .<br/>
To find the neighboring genes of the homologous genes in the databases run this command:<br/>
`python update_neighbor_genes.py -d path -r -test` to test the file for 5 samples or you can directly run:<br/>
`python update_neighbor_genes.py -d path -r -run`. <br/>
It will update/create the `neighbor_genes.json` file in the `processed` directory.<br/>
<br/>**Prepare Synteny Matrices:**
The aim is to get all the alignments of the genes with respect to one another. This step will create the synteny matrix for each record in the datset where the *i,jth* entry of the matrix is the alignment of *ith* gene with *jth* gene. We use both local and global alignments. Also to provide extra features each *ith* gene is reversed and again the alignment is taken again with the *jth* gene. This solves two purposes:<br/>
 1.It makes sure that the `+` strand gene sequence is always aligned with `+` strand gene sequence and vice-versa.<br/>
 2.It provides an extra depth to the synteny matrix which helps in better recognition of features.<br/>

To prepare synteny matrices run this command:<br/>
`python prepare_synteny_matrices.py nos` where:<br/>
`nos` is an integer specifying the number of rows that you want to sample from each file to prepare the dataset. 
This file assumes that:<br/>
`neighbor_genes.json` is present in the `processed` sub-directory.<br/>
homology databases are present in the `data_homology` sub-directory<br/>.
all the `cds` sequence files are present in the `geneseq` directory, but this is not mandatory if you want to directly access all the sequences over the `REST-API`. It accesses all the not found gene sequences from the `REST-API` after reading all the files.<br/>
All the synteny matrices are created and stored in `processed/synteny_matrices` directory in `.npy` format.<br/>
<br/> **Extract Other Features:**<br/>
This will extract other features from the species trees.
Run this command:<br/>
`python prepare_other_features.py`<br/>.
After the completion of this process you should get a binary file named `dataset` in the directory.<br/>

**TILL NOW WE HAVE PREPARED ONLY THE DATASET CONTAINING POSITIVE SAMPLES.**  <br/>
<img src="https://nzspo96lr1-flywheel.netdna-ssl.com/wp-content/uploads/2017/09/tenor.gif" alt="I know, right!!" /> 
<br/>
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;*(I know right)* 

## TO PREPARE NEGATIVE DATASET:
This dataset will contain the gene pairs which are not homologous to each other. <br/>
**Sample the Negative Dataset:**<br/>
Negative dataset takes one gene from the homology database and another gene from the `gtf` files which is not present in any of the homology databases read. This ensures that the two genes selected belong to mutually exclusive sets and hence are not homologous to each other.<br/>
To sample the negative dataset run this command:<br/>
`python prepare_negative_dataset.py nos random_seed` where<br/>
`nos`:number of rows to sample from the dataset.<br/>
`random_seed`:a random number to seed the random number generator.<br/>

**Update Neigbor Genes**<br/>
(Again!!!, ¯\\_(ツ)_/¯). To update the neighbor genes with the new sampled dataset use this command:<br/>
`python update_neighbor_genes.py` <br/>


