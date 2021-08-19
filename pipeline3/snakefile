import numpy as np
import pandas as pd
from datetime import date

the_date = date.today().strftime("%b-%d-%Y")
the_date = "Aug-11-2021" # The name used to produce the output file
# get the list of species names for the homology databases
configfile:
    "config.json"

# sp_path = "database_species_intersection.txt"
species = np.loadtxt(config["species_file"], dtype="str")
pep_paths = pd.Series(np.loadtxt("config/peptide_paths.txt", dtype="str"))
homo_paths = pd.Series(np.loadtxt("config/homology_database_paths.txt", dtype="str"))
gtf_paths =  pd.Series(np.loadtxt("config/gtf_paths.txt", dtype="str"))
models = ["simple_model","ensemble_model","split_feature_model"]
rule all:
	input:
		#expand( config["Diamond_alignments_path"] +"/{SPECIES}.outs.tsv.gz", SPECIES = species)
		# expand("/nfs/production/flicek/ensembl/compara/amarshall/Data_Storage/Longest_Gene_Fasta/{SPECIES}_longest_pep.fa",  SPECIES = species),
		# expand(config["species_neighbour_outdir"] + "/{SPECIES}_all_neighbours.json", SPECIES = species )
		#expand(config["out_dir"] + "/synteny_matrices/{SPECIES}_global2.txt", SPECIES = species[-5:] ),
		# expand(config["out_dir"] + "/Pfam_matrices/{SPECIES}_Pfam.npy", SPECIES = species),
		# expand(config["Longest_pep_fasta_path"] + "/{SPECIES}_longest_pep_rev.fa", SPECIES=species),
		# expand(config["out_dir"] + "/Pfam_matrices/{SPECIES}_negative_Pfam.npy", SPECIES = species),
		# expand(config["Diamond_alignments_path"] + "/{SPECIES}_rev.outs.tsv.gz", SPECIES=species),
		# expand(config["Diamond_alignments_path"] + "/{SPECIES}/{SPECIES}_rev_genes_list.txt", SPECIES=species)
		# expand(config["samples_neighbour_outdir"] + "/{SPECIES}_negative_homology_sample_neighbour.json", SPECIES=species)
		# expand(config["out_dir"] + "/synteny_matrices/{SPECIES}_global2.npy", SPECIES = species),
		# expand(config["out_dir"] + "/negative_synteny_matrices/{SPECIES}_negative_global2.npy", SPECIES = species)
		# expand(config["Diamond_alignments_path"] + "/{SPECIES}/{SPECIES}_genes_list.txt", SPECIES = species)
		# config["Final_data"] + "/" + the_date + "_big_final.csv",
		expand(config["Model_outputs"] + "/{MODEL}/keras_metadata.pb", MODEL = models),
		expand(config["Model_outputs"] + "/model_evaluations" + "/{MODEL}_comparison_matrix.png", MODEL=models)
		# expand(config["Diamond_alignments_path"] + "/{SPECIES}/{SPECIES}_genes_list.txt", SPECIES=species)

# function to map the species to the path of the peptide fasta file
def pep_path_map(wildcards):
	return pep_paths[pep_paths.str.contains(wildcards.SPECIES)].values[0]

rule select_longest:
	input:
		pep_path_map
	output:
		config["Longest_pep_fasta_path"] + "/{SPECIES}_longest_pep.fa"
	threads:
		4
	resources:
		memory = 8000,
		reports_path = config["reports_out_base_path"] + "/select_longest/",
		gpu = ""
	shell:
		"""
		python  scripts/GeneReferencePetpideSelection.py --seq {input} --out_fasta {output}
		"""
rule reverse_fasta:
	input:
		config["Longest_pep_fasta_path"] + "/{SPECIES}_longest_pep.fa"
	output:
		config["Longest_pep_fasta_path"] + "/{SPECIES}_longest_pep_rev.fa"
	threads:
		1
	resources:
		memory = 1000,
		reports_path = config["reports_out_base_path"] + "/reverse_fasta/",
		gpu = ""
	shell:
		"""
		python  scripts/reverse_fasta.py --input_fasta {input} --out_fasta {output}
		"""


rule Diamond_alignment:
	input:
		config["Longest_pep_fasta_path"] + "/{SPECIES}_longest_pep.fa"
	output:
		config["Diamond_alignments_path"] + "/{SPECIES}.outs.tsv.gz"
	params:
		diamond_db_path = config["Diamond_Pfam_db_path"]
	threads:
		16
	resources:
		reports_path = config["reports_out_base_path"] + "/Diamond_alignment/",
		memory = 10000, # LFS memory in MBs
		blocks = 2, # Diamond requires memory ~ blocks * 6Gb
		gpu = ""
	shell:
		"""
		diamond blastp -q {input} -d {params.diamond_db_path} -o {output} --very-sensitive \
		--block-size {resources.blocks} --threads {threads} --verbose --max-target-seqs 250 --compress 1\
		--outfmt 6 qseqid qlen sseqid slen qstart qend sstart send evalue

 		"""

rule rev_Diamond_alignment:
	input:
		config["Longest_pep_fasta_path"] + "/{SPECIES}_longest_pep_rev.fa"
	output:
		config["Diamond_alignments_path"] + "/{SPECIES}_rev.outs.tsv.gz"
	params:
		diamond_db_path = config["Diamond_Pfam_db_path"]
	threads:
		16
	resources:
		reports_path = config["reports_out_base_path"] + "/rev_Diamond_alignment/",
		memory = 10000, # LFS memory in MBs
		blocks = 2,  # Diamond requires memory ~ blocks * 6Gb
		gpu = ""
	shell:
		"""
		diamond blastp -q {input} -d {params.diamond_db_path} -o {output} --very-sensitive \
		--block-size {resources.blocks} --threads {threads} --verbose --max-target-seqs 250 --compress 1\
		--outfmt 6 qseqid qlen sseqid slen qstart qend sstart send evalue

 		"""


def gtf_path_map(wildcards):
	return gtf_paths[gtf_paths.str.contains(wildcards.SPECIES)].values[0]

def neighbour_path_map(wildcards):
	gtf_path = gtf_paths[gtf_paths.str.contains(wildcards.SPECIES)]
	#transform this into the output path of the neighbours
	return config["species_neighbour_outdir"] + gtf_path.str.split("/").str[-1].str.split(".gtf").str[0] + ".json"

rule species_neighbour:
	input:
		gtf_path_map
	output:
		config["species_neighbour_outdir"] + "/{SPECIES}_all_neighbours.json"
	params:
		num_neighbours = config["num_neighbours"],
		out_dir = config["species_neighbour_outdir"]
	threads:
		2
	resources:
		reports_path = config["reports_out_base_path"] + "/species_neighbour/",
		memory = 2000,
		gpu = ""
	shell:
		"""
		mkdir -p {params.out_dir}
		python scripts/geneNeighbourFormater.py --gtf {input} \
		--nb_neighbours {params.num_neighbours} --out_prefix {params.out_dir} \
		--species {wildcards.SPECIES}
		"""

rule split_diamond_alignment:
	input:
		config["Diamond_alignments_path"] + "/{SPECIES}.outs.tsv.gz"
	output:
		config["Diamond_alignments_path"] + "/{SPECIES}/{SPECIES}_genes_list.txt"
	threads:
		10
	resources:
		reports_path = config["reports_out_base_path"] + "/Diamond_Split/",
		memory = 50000,
		gpu = ""
	params:
		out_dir = config["Diamond_alignments_path"] + "/{SPECIES}/"
	shell:
		"""
		bash scripts/split_diamond_alignment.sh {input} {output} {threads} {params.out_dir} {wildcards.SPECIES}
		"""

rule split_rev_diamond_alignment:
	input:
		config["Diamond_alignments_path"] + "/{SPECIES}_rev.outs.tsv.gz"
	output:
		config["Diamond_alignments_path"] + "/{SPECIES}/{SPECIES}_rev_genes_list.txt"
	threads:
		10
	resources:
		reports_path = config["reports_out_base_path"] + "/rev_Diamond_Split/",
		memory = 2000,
		gpu = ""
	params:
		out_dir = config["Diamond_alignments_path"] + "/{SPECIES}/"
	shell:
		"""
		bash scripts/split_rev_diamond_alignment.sh {input} {output} {threads} {params.out_dir} {wildcards.SPECIES}
		"""

def homo_path_map(wildcards):
	return homo_paths[homo_paths.str.contains(wildcards.SPECIES)].values[0]

rule select_samples:
	input:
		config["species_neighbour_outdir"] + "/{SPECIES}_all_neighbours.json",
		homology_database = homo_path_map
	output:
		# The species neighbours
		config["samples_neighbour_outdir"] + "/{SPECIES}_sample_neighbour.json",
		# The species homolog neighbours
		config["samples_neighbour_outdir"] + "/{SPECIES}_homology_sample_neighbour.json",
		# The species neighbours
		config["samples_neighbour_outdir"] + "/{SPECIES}_negative_sample_neighbour.json",
		# The species homolog neighbours
		config["samples_neighbour_outdir"] + "/{SPECIES}_negative_homology_sample_neighbour.json"		
	params:
		out_dir = config["samples_neighbour_outdir"],
		n_samples = config["samples_per_species"],
		n_species = config["num_homolog_species"],
		sp_names = config["sp_names"],
		dist_matrix = config["dist_matrix"],
		neighbour_dir = config["species_neighbour_outdir"]
	threads:
		2
	resources:
		reports_path = config["reports_out_base_path"] + "/select_samples/",
		memory = 2000, # LFS memory in MBs
		# Diamond requires memory ~ blocks * 6Gb
		gpu = ""
	shell:
		"""
		mkdir -p {params.out_dir}

		python scripts/sample_gene.py --homo {input.homology_database} --out {params.out_dir} \
		--species {wildcards.SPECIES} --n_samples {params.n_samples} --n_species {params.n_species} \
		--sp_names {params.sp_names}  --dist_matrix {params.dist_matrix} --neighbours {params.neighbour_dir}
		"""

rule synteny_matrix:
	input:
		# The species neighbours
		config["samples_neighbour_outdir"] + "/{SPECIES}_sample_neighbour.json",
		# The species homolog neighbours
		config["samples_neighbour_outdir"] + "/{SPECIES}_homology_sample_neighbour.json"
	output:
		config["out_dir"] + "/synteny_matrices/{SPECIES}_global1.npy",
		config["out_dir"] + "/synteny_matrices/{SPECIES}_global2.npy",
		config["out_dir"] + "/synteny_matrices/{SPECIES}_local1.npy",
		config["out_dir"] + "/synteny_matrices/{SPECIES}_local2.npy"
	params:
		out_dir = config["out_dir"] + "/synteny_matrices",
		samples_dir = config["samples_neighbour_outdir"],
		longest_fasta_dir = config["Longest_pep_fasta_path"]
	threads:
		2
	resources:
		reports_path = config["reports_out_base_path"] + "/synteny_matrices/",
		memory = 500,
		gpu = ""
	shell:
		"""
		mkdir -p {params.out_dir}
		python scripts/synteny_matrix.py --species {wildcards.SPECIES} --samples_dir {params.samples_dir} \
		--longest_fasta_dir {params.longest_fasta_dir} --out_dir {params.out_dir}
		"""

rule negative_synteny_matrix:
	input:
		# The species neighbours
		config["samples_neighbour_outdir"] + "/{SPECIES}_negative_sample_neighbour.json",
		# The species homolog neighbours
		config["samples_neighbour_outdir"] + "/{SPECIES}_negative_homology_sample_neighbour.json"
	output:
		config["out_dir"] + "/negative_synteny_matrices/{SPECIES}_negative_global1.npy",
		config["out_dir"] + "/negative_synteny_matrices/{SPECIES}_negative_global2.npy",
		config["out_dir"] + "/negative_synteny_matrices/{SPECIES}_negative_local1.npy",
		config["out_dir"] + "/negative_synteny_matrices/{SPECIES}_negative_local2.npy"
	params:
		out_dir = config["out_dir"] + "/negative_synteny_matrices",
		samples_dir = config["samples_neighbour_outdir"],
		longest_fasta_dir = config["Longest_pep_fasta_path"]
	threads:
		2
	resources:
		reports_path = config["reports_out_base_path"] + "/negative_synteny_matrices/",
		memory = 500,
		gpu = ""
	shell:
		"""
		mkdir -p {params.out_dir}
		python scripts/negative_synteny_matrix.py --species {wildcards.SPECIES} --samples_dir {params.samples_dir} \
		--longest_fasta_dir {params.longest_fasta_dir} --out_dir {params.out_dir}
		"""

rule pfam_matrix:
	input:
		# The species neighbours
		config["samples_neighbour_outdir"] + "/{SPECIES}_sample_neighbour.json",
		# The species homolog neighbours
		config["samples_neighbour_outdir"] + "/{SPECIES}_homology_sample_neighbour.json",
		# config["Diamond_alignments_path"] + "/{SPECIES}.outs.tsv.gz",
		config["Diamond_alignments_path"] + "/{SPECIES}/{SPECIES}_genes_list.txt"
	output:
		config["out_dir"] + "/Pfam_matrices/{SPECIES}_Pfam.npy"		
		# config["out_dir"] + "/synteny/" + species + "_local1.txt",
		# config["out_dir"] + "/synteny/" + species + "_local2.txt"
	params:
		out_dir = config["out_dir"] + "/Pfam_matrices",
		samples_dir = config["samples_neighbour_outdir"],
		Pfam_dir = config["Diamond_alignments_path"]
	threads:
		40
	resources:
		reports_path = config["reports_out_base_path"] + "/Pfam_matrices/",
		memory = 4000,
		gpu = ""
	shell:
		"""
		mkdir -p {params.out_dir}
		python scripts/Pfam_matrix_parallell_split.py --species {wildcards.SPECIES} --samples_dir {params.samples_dir} \
		--Pfam_dir {params.Pfam_dir} --out_dir {params.out_dir} --threads {threads}
		"""


rule negative_pfam_matrix:
	input:
		# The species neighbours
		config["samples_neighbour_outdir"] + "/{SPECIES}_negative_sample_neighbour.json",
		# The species homolog neighbours
		config["samples_neighbour_outdir"] + "/{SPECIES}_negative_homology_sample_neighbour.json",
		# config["Diamond_alignments_path"] + "/{SPECIES}.outs.tsv.gz",
		config["Diamond_alignments_path"] + "/{SPECIES}/{SPECIES}_genes_list.txt"
	output:
		config["out_dir"] + "/Pfam_matrices/{SPECIES}_negative_Pfam.npy"		
	params:
		out_dir = config["out_dir"] + "/Pfam_matrices",
		samples_dir = config["samples_neighbour_outdir"],
		Pfam_dir = config["Diamond_alignments_path"]
	threads:
		40
	resources:
		reports_path = config["reports_out_base_path"] + "/negative_Pfam_matrices/",
		memory = 4000,
		gpu = ""
	shell:
		"""
		mkdir -p {params.out_dir}
		python scripts/negative_Pfam_matrix_parallel_split.py --species {wildcards.SPECIES} --samples_dir {params.samples_dir} \
		--Pfam_dir {params.Pfam_dir} --out_dir {params.out_dir} --threads {threads}
		"""

the_date = date.today().strftime("%b-%d-%Y")

rule format_data:
	input:
		expand(config["out_dir"] + "/Pfam_matrices/{SPECIES}_Pfam.npy", SPECIES=species),
		expand(config["out_dir"] + "/Pfam_matrices/{SPECIES}_negative_Pfam.npy", SPECIES=species),
		expand(config["out_dir"] + "/negative_synteny_matrices/{SPECIES}_negative_global1.npy", SPECIES=species),
		expand(config["out_dir"] + "/negative_synteny_matrices/{SPECIES}_negative_global2.npy", SPECIES=species),
		expand(config["out_dir"] + "/synteny_matrices/{SPECIES}_global1.npy", SPECIES=species),
		expand(config["out_dir"] + "/synteny_matrices/{SPECIES}_global2.npy", SPECIES=species)
	output:
		config["Final_data"] + "/" + the_date + "_big_final.npy",
		config["Final_data"] + "/" + the_date + "_big_final.csv"
	params:
		synt = config["out_dir"] + "/synteny_matrices",
		neg_synt = config["out_dir"] + "/negative_synteny_matrices",
		out_dir = config["out_dir"],
		samples_dir = config["samples_neighbour_outdir"],
		Pfam_dir = config["out_dir"] + "/Pfam_matrices",
		date = the_date,
		Final_data = config["Final_data"]
	threads:
		2
	resources:
		reports_path = config["reports_out_base_path"] + "/combine_data/",
		memory = 80000,
		gpu = ""
	shell:
		"""
		mkdir -p {params.out_dir}
		python scripts/combine_data.py --samples_dir {params.samples_dir} --Pfam_dir {params.Pfam_dir} \
		--synt {params.synt} --neg_synt {params.neg_synt} --out_dir {params.out_dir} --date {params.date} --Final_data {params.Final_data}
		"""

rule train_model: 
	input:
		config["Final_data"] + "/" + the_date + "_big_final.npy",
		config["Final_data"] + "/" + the_date + "_big_final.csv"
	output:
		config["Model_outputs"] + "/{MODEL}/{MODEL}_accuracy_curves.png",
		config["Model_outputs"] + "/{MODEL}/{MODEL}_loss_curves.png",
		config["Model_outputs"] + "/{MODEL}/keras_metadata.pb"
	params:
		out_dir = config["Model_outputs"]
	threads:
			4
	resources:
		reports_path = config["reports_out_base_path"] + "/Train_models/",
		memory = 10000,
		# gpu = ""
		gpu = '-q gpu -gpu "num=1:j_exclusive=no"'

	shell:
		"""
		mkdir -p Model_outputs
		python scripts/Deep_Learning/on_disc_trainer.py --npy {input[0]} --csv {input[1]} \
		--model_name {wildcards.MODEL} --out_dir {params.out_dir}
		"""

rule evaluate_model:
	input:
		config["Final_data"] + "/" + the_date + "_big_final.npy",
		config["Final_data"] + "/" + the_date + "_big_final.csv",
		config["Model_outputs"] + "/{MODEL}/keras_metadata.pb" # check for the trained model
	output:
		config["Model_outputs"] + "/model_evaluations" + "/{MODEL}.png",
		config["Model_outputs"] + "/model_evaluations" + "/{MODEL}_comparison_matrix.png"
	params:
		out_dir = config["Model_outputs"] + "/model_evaluations",
		models_dir = config["Model_outputs"]
	threads:
			4
	resources:
		reports_path = config["reports_out_base_path"] + "/Evaluate_models/",
		memory = 5000,
		gpu = ""
		# gpu = '-q gpu -gpu "num=1:j_exclusive=no"'

	shell:
		"""
		mkdir -p {params.out_dir}
		python scripts/Deep_Learning/evaluator.py --npy {input[0]} --csv {input[1]} \
		--model_name {wildcards.MODEL} --out_dir {params.out_dir} --models_dir {params.models_dir}
		"""