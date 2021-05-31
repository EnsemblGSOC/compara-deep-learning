# Project Notebook

This is markdown file is used as a notebook of the project to log what has/hasn't been done, document problems, and track progress

Mon May 31 06:33:39 BST 2021

Encountered problem with the neighbour genes. I think this is because I didn't select a large enough number of genes from each homology database. I chose 2 originally. There are 20,000 genes so it isn't topo surprising this didn't work. I am upping this to 1,000 and hoping this works.

I Fixed this and was able to produce the relevant files. I discovered the above was not the main source of the problem however. What I need to do is run hmmerscan, which is poorly described in the current github description. I will update this in due course. I have installed hmmer using the mamba approach here and found the manual. I believe the correct next task is to follow the tutorial in their online manual. To get the hmmerscan function to run correctly you need hmmdb, which I think means run a database. I think the correct database is the pfam one, but right now this is just a hunch. I will continue from here tomorrow.




