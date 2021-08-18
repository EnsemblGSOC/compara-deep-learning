pipeline_dir="/hps/software/users/ensembl/repositories/compara/amarshall/compara-deep-learning/pipeline3/"
echo $pipeline_dir

cd $pipeline_dir
export OMP NUM THREADS=16
snakemake --dag | dot -Tpng > dag.png
snakemake -s snakefile --use-conda --cluster \
'bsub -e {resources.reports_path} -o {resources.reports_path} -n {threads} -M {resources.memory} {resources.gpu}' \
--jobs 220 --latency-wait 30