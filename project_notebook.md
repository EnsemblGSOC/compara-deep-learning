# Project Notebook

This is markdown file is used as a notebook of the project to log what has/hasn't been done, document problems, and track progress

Mon May 31 06:33:39 BST 2021

Encountered problem with the neighbour genes. I think this is because I didn't select a large enough number of genes from each homology database. I chose 2 originally. There are 20,000 genes so it isn't too surprising this didn't work. I am upping this to 1,000 and hoping this works.

I Fixed this and was able to produce the relevant files. I discovered the above was not the main source of the problem however. What I need to do is run hmmerscan, which is poorly described in the current github description. I will update this in due course. I have installed hmmer using the mamba approach here and found the manual. I believe the correct next task is to follow the tutorial in their online manual. To get the hmmerscan function to run correctly you need hmmdb, which I think means run a database. I think the correct database is the pfam one, but right now this is just a hunch. I will continue from here tomorrow.

Tue Jun  1 05:44:21 BST 2021

The pfam ftp has reference databases to perform the alignments against: http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/.
From this, I am downloading Pfam-A.seed.gz. I believe this is the PFAM database against which protein sequences are aligned. I am running the command hmmbuild Pfam-A.hmm Pfam-A.seed to build the Pfam profile to perform sequence alignments. I then ran hmmpress to compress this, as stated in the hmmer manual. Finally, we can run hmmscan. (And it all actually runs. AMAZING!)

Reading through the code for the cmap, it looks like it's list of dictionaries. Maybe this data structure could be improved by using arrays and the appropriate indexing instead. It depends on how much of a bottle-neck this really turns out to be in the end.

Wed Jun  2 08:04:06 BST 2021

I am running hmmrscan on the positive file. However, this is really really slow as stated in the project description. I am going to run it on the positive examples, however, I am going to download GSoC 2019's domtblout file from http://ftp.ebi.ac.uk/pub/databases/ensembl/mateus/gsoc_2019/pfam_scan/pro_seq_negative/pro_seq_negative.fa.domtblout.gz. The negative examples that I run has far more examples than the positive version I am doing myself and would take probably days to align so downloading is expedient. I will run a full pipeline from start to end at a future date.

Thu Jun  3 04:42:57 BST 2021
I Building on what I said yesterday, I believe that I have managed to successfully run the entire pre-processing pipeline for the 50k negative examples data-set. Now I need to understand the NN structure and understand the different inputs to the traning script. 

Mon Jun 14 14:50:13 BST 2021
I re-wrote the pipeline to work with the the full species list from all of the different databases. I then re-ran most but not all of the pipeline.

Having spoken to David and Jorge, I should attempt to write a database of all pfam alignments to disc. There are two approaches to this. The first is to figure out how to get DIAMOND to play ball with the pfam hmm databases. The second is to generate this as a one time only thing using hmmer, version 3. I need to ask the team about what settings of pfam they really care about being used. To do this I should compare this to the outputs that were obtained from HMMER in the previous project and what we can get from DIAMOND instead of this. Only after I have done this should I generate the full list of Alignments for every species out there.

Tue Jun 22 05:08:59 BST 2021

Last week was an underperforming week due to the severely poor performance of the Imperial hpc and my inability to flexibly adapt to that situation. The full Pfam-A.fasta file takes far too long to perform alignments on my local machine. I am therefore exploring methods to subsample this fasta file to have a prototype pipeline with diamond up and running. 

To subsample 10,000 random sequences from the full Pfam-A.fasta file I am using the tool seqtk with the command: 

seqtk sample -s seed=42 Pfam-A.fasta.gz 10000 | gzip > Pfam-A_subsample.fasta.gz

(For more info on seqtk see here https://github.com/lh3/seqtk. I installed using Mamba)

Thu Jun 24 07:45:14 BST 2021
I have added modified the Pfam overlap scoring function to a Jaccard index. This is likely sub-optimal but allows for the project to move forwards in the meantime. The reason for the change is that if gene aligned to a Pfam domain in two separate locations the old python code threw an error because it didn't anticipate multiple alignments. The Jaccard is a simplistic solution to that problem. Also, the old scoring function used the gene with the most Pfam domains as the denominator for the score. This seems asymmetric and weird. Jaccard doesn't suffer this problem.

This implementation has at least successfully run.

Fri Jun 25 10:04:12 BST 2021

It looks like the positive examples were processed to make Pfam Jaccard matrcies successfully. I am just waiting on the same for the negative examples list which is much much larger. 

I have just realised the Finalize data script doesn't actually utilist the PFAM matrices in anyway. This needs modifcation!!!

I have sketech the outlines of a conv net in a ipynb. This was based on the old HMMER alingment data set that I made. It achieved high accuracy due to the severe class imbalances caused by how large the negative dataset is. Basically, just be on a negative example and you would do well. 

Sat Jun 26 07:35:23 BST 2021

I have successfully been able to get the PFAM matrices into the finalised dataset using the diamond alignments from the data. I am however confused by the end stage of the process. My impression was that PFAM matrices from diamond, and synteny matrices were to be used to encode some form of relation ship and that this was to be used to form the nn.

I have begun to iterate on simple versions of conv nets which take in the PFAM matrices as well as the Synteny matrices. The positive examples dataset is however limited in it's size Inhibiting performance.

I will now perform an upscale of the pipeline to produce a bigger training set. We will multiply the number of training examples by around 50, and scale up the Diamond database to have around 1 million PFAM domains instead of the 10,000 I did before. As always, the random seed was 42.

Wed Jul  7 15:19:29 BST 2021

Significant progress has been made! The EBI hpc is FANTASTIC. I have migrated a good chunk of the workflow over to their system. It runs very smoothly. I have started to do a re-write of the pipeline so that it isn't horifically slow. This includes:
- performing Diamond alignments for the full Pfam-A database for every protein in a species. Specifically the longest version of every protein.
- Utilising David's far faster scripts for extracting neighbours from the pipeline.
- A new and faster sampling approach that makes use of pandas.

The LSF hpc system works flawlessly for job submissions with snakemake and the integration was one of the nicest experiences I've had with anything from bioinformatics.

To do:

I have yet to incorporate the pairwise aligner for the synteny matrices. I think David's scripts for this look far more sensible than the previous approaches. I'll be looking to integrate this next. I will be using the Jaccard matric for the Pfam matches between genes just like last time.
