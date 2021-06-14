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
