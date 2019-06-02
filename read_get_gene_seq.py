import json
from Bio import SeqIO
import pandas as pd
import os
import gzip

def read_from_multiple_lsy(lsyfl):
    lsy={}
    for f in lsyfl:
        d={}
        with open(f,"r") as file:
            d=dict(json.load(file))
            for t in d:
                lsy[t]=d[t]
    return lsy

#this function updates the given dictionary with the given keys and values list
def create_dict(keys,values,dictionary):
    for i in range(len(keys)):
        if keys[i] not in dictionary:
            dictionary[keys[i]]=values[i]

    return dictionary

#this function maps all the genes to their respective species.
#(Function: when finding the species of any gene we do not need to search the entire dataframe)
def group_seq_by_species(df):
    g_to_sp={}
    sph=list(df.homology_species)
    ghsp=list(df.homology_gene_stable_id)
    sp=list(df.species)
    gsp=list(df.gene_stable_id)
    create_dict(gsp,sp,g_to_sp)
    create_dict(ghsp,sph,g_to_sp)
    return g_to_sp

#this function returns the gene-id and gene-biotype from the description in the fasta file record.
def description_cleaner(description):
    description=description.split()
    t=""
    gbt=""
    for x in description:
        try:
            x=x.split(":")
            if x[0]=="gene":
                t=x[1].split(".")[0]
            if x[0]=="gene_biotype":
                gbt=x[1]
        except:
            return "aa","aa"
    return t,gbt

def read_gene_seq(dirname,s,genes_by_species):
    lof=os.listdir(dirname)#list all the files in the sequences directory
    ftr=[]
    for f in lof:
        if f.split(".")[0] in s:#check whether the species is present in the species to read list. Will skip those species which are not present in the dataframe
            ftr.append(f)
    data={}
    for f in ftr:
        species=f.split(".")[0].lower()
        with gzip.open(dirname+"/"+f,"rt") as file:
            record=SeqIO.parse(file,"fasta")
            for r in record:
                gid,gbt=description_cleaner(r.description)
                if str(gid) not in data and str(gid) in genes_by_species[species] and gbt=="protein_coding":
                    data[gid]=str(r.seq)
    return data

def read_gene_sequences(df,lsy,data_dir,fname):

    """The basic idea here is to create a list/dictionary of all the genes by their species.
       Once the mapping is done, all the respective fasta sequence files are read by Species
       and the CDNA sequences for each gene in the species record are read and stored.
       Thus we don't have to read the same file multiple times."""

    grouped_genes=group_seq_by_species(df)
    gene_by_species_dict={}
    for i in df.homology_species.unique():
        gene_by_species_dict[i]=[]
    for x in lsy:
        species=grouped_genes[x]#get the species
        if x not in gene_by_species_dict[species]:#check if the gene already exists in the species dict or not.
            gene_by_species_dict[species].append(x)
        xl=lsy[x]['b']
        xr=lsy[x]['f']
        for gxl in xl:
            if gxl=="NULL_GENE":
                break
            if gxl not in gene_by_species_dict[species]:
                gene_by_species_dict[species].append(gxl)
        for gxr in xr:
            if gxr=="NULL_GENE":
                break
            if gxr not in gene_by_species_dict[species]:
                gene_by_species_dict[species].append(gxr)

        s=[x for x in gene_by_species_dict if len(gene_by_species_dict[x])!=0]#select those species only whose gene sequences we have to read.
        s=[x.capitalize() for x in s]

        data=read_gene_seq(data_dir,s,gene_by_species_dict)

        with open("processed/"+fname+".json","w") as file:#save the data
            json.dump(data,file)

        return data
