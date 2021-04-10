The aim of this project is to use **Deep Neural-Nets** to predict homology type between a given pair of genes. 

## Requirements:
1. A machine with atleast **8GB of RAM** (although **16-32GB** is recommended. It depends on the no. of homology databases that you are willing to use in the preparation of the dataset), a graphic card for training the deep neural nets. A single GPU machine would suffice. The model can be trained on CPU as well but will be a lot faster if trained on a GPU.<br/>
2. A stable Internet Connection.<br/>
3. A native/virtual python environment. Install the dependencies from `requirements.txt` using :<br/>
 `pip install -r requirements.txt`<br/>
 
## Data Preparation:
The model uses a synteny matrix and some other factors derived from the species tree to make predictions. <br/>
**Download The Required Files:**

In order to download all the needed files run: `python ftpg.py`<br/>

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

All the databases have the same name so you have to rename the files with their respective species names.
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
Negative samples are a non-homologous pair of genes. You can get the negative samples from [here](ftp://ftp.ebi.ac.uk/pub/databases/ensembl/mateus/gsoc_2019/).Download one of your choice :).<br/>
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
