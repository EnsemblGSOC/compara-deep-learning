#!/bin/sh

echo $LSB_JOBINDEX
species=$( head -$LSB_JOBINDEX config/valid_species.txt | tail -1)
echo "running script!"
python scripts/synteny_matrix.py --species $species \
--samples_dir /nfs/production/flicek/ensembl/compara/amarshall/Data_Storage/Sample_Neighbours \
--longest_fasta_dir /nfs/production/flicek/ensembl/compara/amarshall/Data_Storage/Longest_Gene_Fasta \
--out_dir /hps/software/users/ensembl/repositories/compara/amarshall/compara-deep-learning/pipeline3/outs/test

exit 0