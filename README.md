# New Readme

The aim of this project is to use **Deep Neural-Nets** to predict homology type between the give pair of genes. 
This file provides the instructions to replicate the results of this project from scratch.

## Requirements:
1. A machine with atleast **8GB of RAM** (although **16-32GB** is recommended. It depends on the no. of homology databases that you are willing to use in the preparation of the dataset), a graphic card for training the deep neural nets. A single GPU machine would suffice. The model can be trained on CPU as well but the training process will be faster if trained on a GPU.<br/>
2. A stable Internet Connection.<br/>
3. A native/virtual python environment. Install the dependencies from `requirements.txt` using :
 `pip install -r requirements.txt`<br/>

## Step 1: Data Preparation:
The model uses a synteny matrix and some other factors derived from the species tree to make predictions. <br/>
**Download The Required Files:**
So as to prepare data we will need the following files:<br/>
1.All the `gtf` files so as to get the start and end locations of the genes and find their neighboring genes. This is used to create the synteny matrix which helps to see the conserved synteny among the genes.<br/>
2. All the `cds` files in `FAST-A` format. They are required but are not mandatory, if the files are not provided the sequences are directly accessed from the `REST API` but the process can be slow :( . It's better to have all the `cds` files.<br/>
3. Homology Databases of Your Choice. All the databases have the same name so its better to change the name of the files with their respective speicies names. One format that works best is `species_name.tsv.gz`.<br/>

To download the `gtf` and `cds` files `ftpg.py` can be used. This scripts writes the links of all the required `gtf` and `cds` files to `gtf-link.txt` and `seq_link.txt`.<br/> You can use your own script to download the files or just enter `y` when prompted for permission to download the files.<br/> It will automatically download all the files and store them in designated folder. You can manually download each file by copying each link from the files and pasting it in the browser.<br/>

The designated folders to store the files are as follows:<br/>
`data` all the `gtf` files.<br/>
`data_homology` all the homology databases that you want to use.<br/>
`geneseq`To store the `cds` sequence files in `fasta` format.<br/>
You can assign different folder names if you want but it's better to stick to them. <br/>
<br/>**Create Genome Maps:**<br/>
The purpose is to create maps of all the genes present in the `gtf` files with respect to their chromosomes, a map of all the genes belonging to the same chromosome in the given species, a map of all the genes in the given species, a map of all the species whose data has been successfully read.<br/>
To create genome maps run this command:<br/>
`python create_genome_maps.py -d path -r` where:<br/>
`path`: link to the directory where all the `gtf` files exist. If you have used `ftpg.py` then the path is `data`. <br/>
Note: Genome Maps can be downloaded from this [link](https://drive.google.com/open?id=1GjV6dT-Hpf2LWQ-vSpekqqQ7RF_tH8So).<br/>
<br/>**Create/Update Neighbor Genes File:**<br/>
This project uses the measure of conserved synteny to predict the homology type. Therefore, to predict the homology type we need to find the neighboring genes of the given homologous pair of genes. <br/>
This file finds the neighboring genes of all the rows in the given databases and writes it to a file called `processed/neighbor_genes.json` .<br/>
To find the neighboring genes of the homologous genes in the databases run this command:<br/>
`python update_neighbor_genes.py -d path -r -test` to test the file for 5 samples or you can directly run:<br/>
`python update_neighbor_genes.py -d path -r -run`. <br/>
It will update/create the `neighbor_genes.json` file in the `processed` directory.<br/>
