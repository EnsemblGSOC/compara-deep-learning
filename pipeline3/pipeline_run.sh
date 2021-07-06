pipeline_dir="/hps/software/users/ensembl/repositories/compara/amarshall/compara-deep-learning/pipeline3"
echo $pipeline_dir

export OMP NUM THREADS=6
bsub -n 6 <your code>

snakemake --dag | dot -Tpng > dag.png
snakemake -s snakefile --use-conda --cluster 'bsub -e $pipeline_dir -o $pipeline_dir -n {threads} -M {resources.memory}' --jobs 10