The aim of this project is to use **Deep Neural-Nets** to predict homology type between a given pair of genes. 

## Requirements:
1. A machine with atleast **8GB of RAM** (although **16-32GB** is recommended. It depends on the no. of homology databases that you are willing to use in the preparation of the dataset), a graphic card for training the deep neural nets. A single GPU machine would suffice. The model can be trained on CPU as well but will be a lot faster if trained on a GPU.<br/>
2. A stable Internet Connection.<br/>
3. A native/virtual python environment. Install the dependencies from `requirements.txt` using :<br/>
 `pip install -r requirements.txt`<br/>
 
## Data Preparation:
The model uses a synteny matrix and some other factors derived from the species tree to make predictions. <br/>
**Download The Required Files:**

### Reference sequence files

3 file formats are required as pre-requisites for this project: `gtf`, `cds` and `pep` files. These files are

1) Gene transfer format (gtf): Gives details of start & end positions of all genes within a species aswell as gene names  http://m.ensembl.org/info/website/upload/gff.html
2) CDS: gives relative position of basepair in *coding sequence*, ie the bit that goes onto determine the amino acid sequence
3) PEP: Contains the amino-acid strings for each gene

The list of every species available is queried using the ftpg.py script. This script goes to the ensemble ftp and yanks the links to all the `gtf`, `cds` and `pep` files available, and writes those links to the text files in the directory download_links/. There will be one text file for each of the 3 file types.

These files will be used later to ...... *(fill this in)*

The bash script, downloader.sh uses wget to download files, and can be run as follows:

downloader.sh download_links.txt out_dir_path out_log_path

You just need to supply it the list of links to download, the output directory to download the files to and a path to a log file for the download.

The snakefile calls this script 3 times, with the relevant inputs to download the set of `gtf`, `cds` and `pep` files to the directory downloads/

### Homology databases

To run you will also need a set of homology databases.

ftp://ftp.ensembl.org/pub/current_tsv/ensembl-compara/homologies/

If you want to eye ball what species are available here then run:

wget ftp://ftp.ensembl.org/pub/current_tsv/ensembl-compara/homologies/ 

This will download an html file, likely called index.html. If you open this html in a regular web browser you can see the list of species on offer.

If you inspect one of these files, eg  homology_databases/homo_sapiens.homologies.tsv.gz you'll see that the left hand column displays a human gene, with an unfriendly name like ENSG00000271254 . Further along is a column called homology_gene_stable_id which gives the id of a gene like ENSNLEG00000033214. If you go this page --- https://asia.ensembl.org/Multi/Search/ --- you can look up what these genes are. You'll see that ENSNLEG00000033214 is *novel gene (Northern white-cheeked gibbon Gene)*, whilst ENSG00000271254 is a human gene. This matches what is written in the file for the species. The homology tsv file will tell you the type of homology relationship between these species, eg "ortholog_many2many". What each of these types means is detailed here: https://m.ensembl.org/info/genome/compara/homology_types.html. The file also contains a GOC score as well which will further be used as a feature prediction for the NN *(I think this is true, double check)*.


## The Genome Map

So earlier on we downloaded a set of gtf files from many many species, and each one outlines the position of all the genes on the genome of that particular species. This point in the pipeline uses a bunch of scripts, located in the create_genome_maps/ directory, to create a map of all the genes. It's a pickle file so you can read it using pandas.read_pickle. This map contains in it's current instantiation from the previous GSoC (on date Wed Jun  9 06:16:40 BST 2021) is structured as follows:

There is a top level dictionary with the following keys: ['cmap', 'cimap', 'ld', 'ldg', 'a', 'd'].
Behind each of these keys is a list of the same length as the number of species in set of gtf files that you downloaded earlier. So within each of the following dictionaries you'll get one of these per species. Each contains the following:

- genome_maps["cmap"] is a list of dictionaries mapping each gene to each chromosome for each species
- genome_maps["cimap"] (*I am unsure*) gives the list of gene indexes on each chromosome. See genome_maps["ldg"] below. 
- genome_maps["ld"] is a dictionary containing lists of genes for each species
- genome_maps["ldg"] (*I am unsure*) I think it's a dictionary mapping each gene to an arbitrary index  (is a list of dictionaries containing the start positions for each gene on each chromosome)
- genome_maps["a"] I think is a dataframe from version of the gtf file including only the protein coding regions. The indexes match those in genome_maps["ld"]
- genome_maps["d"] maps each species to it's index position in all of the other lists. For example, if you type 

'''
genome_maps["d"]["Homo_sapiens"]
'''
you'll get back

72

If you then go and look at genome_maps["cmap"][72], you'll see the list of all genes and their chromosomes for homosapiens.


The geneome map contains a lot of redundant information at the moment.


Roughly, this map is created by going through every downloaded gtf file as a pandas dataframe and retaining only the genes from the gtf (as opposed to transcripts, or some other genomic feature). We get a list of all these gene dataframes, one for each species, and a dictionary which takes a species name and tells you it's position in the list of gtf dataframes just created (this is done by the read_data.py script's read_data_genome function).  The get_data.py script's get_data_genome function uses this dataframe list and dictionary to...... *(This get's fairly convoluted without any available doc-strings. Complete this description of the cmap!!!)*.


GOC scores
If you look inside of the homology databases you'll find that there's a column called goc score which have values of 0,25,50,100. This is related to the number of different orthologs that nearby neighbouring genes have. See here for more details about its computation  https://asia.ensembl.org/info/genome/compara/Ortholog_qc_manual.html.






Tasks to incorporate. 
- Finish the cmap description
- Understand what pfam is in depth
- Understand loosely what hmmerscan is doing. Consider how best to substitute diamond in hmmerscan's place




It will download the `gtf`, `cds` and `pep` files.

This scripts writes the links of all the required `gtf` and `cds` files to `gtf-link.txt` and `seq_link.txt`. You can use your own script to download the files or just enter `y` when prompted for permission to download the files. It will automatically download all the files and store them in designated folder. You can manually download each files by pasting link from the files in the browser.<br/>

So as to prepare data we will need the following files:<br/>
1.All the `gtf` files to get the start and end locations of the genes and find their neighboring genes. This is used to create the synteny matrix which helps to see the conserved synteny among the genes.<br/>
2. All the `cds` files in `FAST-A` format. They are required but are not mandatory, if the files are not provided the sequences are directly accessed from the `REST API` but the process can be slow :( . It's better to have all the `cds` files.<br/>
3. All the `pep` files in `FAST-A` format. They are required to get the protein sequences to run the pfam scan on.<br/>

**Homology Databases:**<br/>
Additionally, you will need to download the homologies of your choice.
They can be found at: ftp://ftp.ensembl.org/pub/current_tsv/ensembl-compara/homologies/<br/>

e.g.: ftp://ftp.ensembl.org/pub/current_tsv/ensembl-compara/homologies/homo_sapiens/Compara.97.protein_default.homologies.tsv.gz

All the databases have the same name so you have to rename the files with their respective speicies names.
From `Compara.97.protein_default.homologies.tsv.gz` to `species_name.tsv.gz`

**Homology Databases:**<br/>
Additionally, you will need to download the homologies of your choice.
They can be found at: ftp://ftp.ensembl.org/pub/current_tsv/ensembl-compara/homologies/<br/>

e.g.: ftp://ftp.ensembl.org/pub/current_tsv/ensembl-compara/homologies/homo_sapiens/Compara.97.protein_default.homologies.tsv.gz

All the databases have the same name so you have to rename the files with their respective speicies names.
From `Compara.97.protein_default.homologies.tsv.gz` to `species_name.tsv.gz`
All homology files must go into a directory called: `data_homology/`

**Create Genome Maps:**<br/>
The purpose is to create maps of all the genes present in the `gtf` files with respect to their chromosomes, a map of all the genes belonging to the same chromosome in the given species, a map of all the genes in the given species, a map of all the species whose data has been successfully read.<br/>
To create genome maps run this command:<br/>
`python create_genome_maps.py`:<br/>
Note: A precomputed Genome Maps can be downloaded from this [link](https://drive.google.com/open?id=1GjV6dT-Hpf2LWQ-vSpekqqQ7RF_tH8So).<br/>

**Select the Records from each homology database:**<br/>
This step will select the specified no. of records from each of the homology databases on the basis of distant species,GOC score,homology type etc.<br/>
To select the data run:<br/>
`python select_data.py number_of_records_to_be_selected_from_each_file`.<br/>

**Find the Neighbor Genes of the Selected Records:**<br/>
This step will find the neighbor genes of all the selected records from the homology databases and write it to the `processed` directory.<br/>
Run:<br/>
`python neighbor_genes.py`<br/>

**Create Synteny Matrix and Write Protein Sequences:**<br/>
This step will create the synteny matrices and write the protein sequences on a `FAST-A` file named `protein_seq_positive.fa`.You will have to run the `hmmer-scan` on this file to get the `PFAm` domains.<br/>
To create the synteny matrices:<br/>
`python prepare_synteny_matrix.py number_of_threads`.<br/>
Note:The script uses multiprocessing to prepare synteny matrices. Since each thread has it own copy of all the resources its better to run with higher number of threads on a high-RAM machine.

**Process The Negative Dataset:**<br/>
Negative samples are a non-homologous pair of genes. You can get the negative samples from here [link](ftp://ftp.ebi.ac.uk/pub/databases/ensembl/mateus/gsoc_2019/).Download one of your choice :).<br/>
Process the negative set by using:<br/>
`python process_negative.py negative_database_file_name number_of_threads`<br\>

**Run the HMMER scan on the Protein Sequences:**<br/>
The idea is that homologous genes will have overlapping domains. Run the hmmer scan on the `protein_seq_positve.fa` and `pro_seq_negative.fa` with the `-domtblout` option.

**Parse the PFAm domain files:**<br/>
This file will parse the hmmer scan database. You will have to parse both the positive samples and the negative sample database.<br/>
To parse, run:<br/>
`python pfam_parser.py domtblout_file_name_positive domtblout_file_name_negative`

**Create PFAM matrices:**<br/>
This step will create the pfam matrices. This might take some time...<br/>
Run:<br/>
`python pfam_matrix.py negative_samples_50K.txt`<br/>

**Finalize the Dataset:**<br/>
This step combines everything and finalizes the dataset by reading the processed factors and extracting some basic features from the species tree. 
Run:<br/>
`python finalize_dataset.py name_of_negative_database_you_earlier_processed`<br/>

IF YOU DID EVERYTHING RIGHT YOU SHOULD SEE A FILE NAMED `dataset` IN THE SAME DIRECTORY.<br/>

## Train the Model:
This is where it gets interesting. You are gonna train your own model architecture or you can use the one given in the `model.py` script. You can directly change the model as well in the `model.py`.<br/>
To train the model, run:<br/>
`python train.py model_name negative_start_composition negative_end_composition epochs learning_rate learning_rate_decay no_of_samples batch_size`<br/>

## Predictions:
To make predictions you need to have the prediction files in a pre-defined format like this [file](ftp://ftp.ebi.ac.uk/pub/databases/ensembl/mateus/gsoc_2019/balanced_random_mix_ortholog_paralog_negative.txt.gz). All the fields have to tab seperated and in the same order.<br/>

**Get the prediction files:**<br/>
Create a new directory with any name of your choice in the code directory and paste all the files on which you want to make the predictions inside it.<br/>
Run:<br/>
`python pfam_folder_pred.py directory_name`<br/>
This will read all the files on the directory and write all the protein sequences on a FAST-A file named `prediction_directory_name.fa`.
Run hmmer scan on that file.<br/>

**Parse the PFAM file**<br/>
This file parses the PFAM file and creates some maps.<br/>
Run:<br/>
`python pfam_db_parser.py domtblout_file_name`<br/>

**Now you have all the resources required to Make predictions on the required file.**<br/>
Copy the file you want to make predictions on from the directory that you created earlier to the code directory and Run:<br/>
`python prediction_pfam.py file_name_with_extension model_name no_of_models number_of_threads start end name domtblout_file_name weights`.<br/>
where:<br/>
`start end`:the locations from which you want to make predictions in the file.<br/>
`name`:Since you can run multiple predcitions at the same time this serves as a unique identifier for the temporary files being created.<br/>
`domtblout_file_name`:It is the name of the file that you get after running the hmmer scan on the fast-a files.<br/>
`weights`:`e` if you want to give equal weight to all the models, `w` followed by the same number of float literals as number of models parameter to assign specific weights. 
