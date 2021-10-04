from typing import List, Dict
import sys
import argparse
import json
import os
from gtf_adaptor import GTFAdaptor


def fetchNeighbourDictionary(fileName: str, nbNghbr: int) -> Dict:
    """
    """

    gtf_adaptor = GTFAdaptor(fileName)
    gtf_adaptor.load_genes("protein_coding")
    results = {}
    for (geneId, val) in gtf_adaptor.genes.items():
        up_genes, down_genes = gtf_adaptor.getNeighbourGenes(geneId, nbNghbr)
        up_gene_id = []
        for gene in up_genes:
            if gene is None:
                up_gene_id.append("NA")
            else:
                up_gene_id.append(gene["stable_id"])
        down_gene_id = []
        for gene in down_genes:
            if gene is None:
                down_gene_id.append("NA")
            else:
                down_gene_id.append(gene["stable_id"])
        results[geneId] = [up_gene_id, down_gene_id]
        # results[geneId]["upstream"] = up_gene_id
        # results[geneId]["downstream"] = down_gene_id
    return results


def listGtfFiles(directory: str) -> List:
    gtf_files = []
    fileNames = os.listdir(directory)
    for fileName in fileNames:
        if fileName[-4:] == ".gtf":
            gtf_files.append(directory + "/" + fileName)
    return gtf_files


def main(args):
    nbNghbr = args.nb_neighbours
    outPrefix = args.out_prefix
    gtfFiles = []
    gtfFiles.append(args.gtf)

    # if a directory  then get the list of gtf file
    if args.dir != "":
        gtfFiles = listGtfFiles(args.dir)

    # fetch neighbor in all gtf of the list
    for gtfFiles in gtfFiles:
        print("fetch neighbours in " + gtfFiles)
        neighbours = fetchNeighbourDictionary(gtfFiles, nbNghbr)
        neighbourFile = (
            outPrefix + "/" + args.species + "_all_neighbours.json"
        )
        with open(neighbourFile, "w") as file_handler:
            file_handler.write(json.dumps(neighbours))


########################################################################################################################

parser = argparse.ArgumentParser()
parser.add_argument("--gtf", help="gtf annotation file", default="")
parser.add_argument("--dir", help="directory all gtf annotation file", default="")
parser.add_argument(
    "--nb_neighbours", type=int, help="nb of neighbour gene to bve retrieved"
)
parser.add_argument(
    "--species", type=str, help="Species name of homology db being processed"
)
parser.add_argument(
    "--out_prefix", help="outfile with the neighbour information", default="."
)
args = parser.parse_args()

if args.gtf == "" and args.dir == "":
    print("--gtf or --dir need to be set. Currently none of them are set")
    sys.exit()
if args.gtf != "" and args.dir != "":
    print(
        " --gtf and --dir are both exclusive so you cannot have them set on the same time"
    )
    sys.exit()

main(args)

