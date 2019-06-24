import os
import pandas as pd
import gzip
import sys
import progressbar

def clear_data(x):
    x=x.split()
    try:
        x=x[1]
    except:
        c=0
        #print(x)
    x=x[1:-1]
    return x

def read_data_genome(dir_name,a,dict_ind_genome):
    lf=os.listdir(dir_name)
    if len(lf)==0:
        print("No files in the data directory!!!!!!")
        sys.exit(1)
    colname=["Chr","source","feature","start","end","score","strand","frame","attribute"]
    print("Going to read data:")
    for x in progressbar.progressbar(range(len(lf))):
        data_gene=pd.read_csv(dir_name+"/"+lf[x],compression='gzip',sep='\t',comment='#',header=None,names=colname)
        #print(data_gene.head)
        data_gene=data_gene[data_gene["feature"]=="gene"]
        tmp=data_gene["attribute"].str.split(";",expand=True)
        tmp=tmp.iloc[:,:5]
        data_gene[["gene_id","gene_version","gene_name","gene_source","gene_biotype"]]=tmp
        data_gene=data_gene.drop("attribute",axis=1)
        #print(data_gene[0:10])
        try:
            for y in ["gene_version","gene_name","gene_source","gene_biotype","gene_id"]:
                data_gene[y]=data_gene[y].apply(clear_data)
        except:
            continue
        #print(data_gene[0:10])        
        data_gene=data_gene[(data_gene['gene_biotype']=='protein_coding') | (data_gene['gene_source']=='protein_coding')]
        #print(data_gene[data_gene["gene_id"]=="ENSPMGG00000022088"])
        a.append(data_gene)
        n=lf[x].split(".")[0]
        dict_ind_genome[n]=len(a)-1
    return a,dict_ind_genome

def read_data_homology(dir):
    a_h=[]
    d_h={}
    lf=os.listdir(dir)
    if len(lf)==0:
        print("No Files in the Directory!!!!!!!")
        sys.exit(1)
    for x in lf:
        data=pd.read_csv(dir+"/"+x,compression='gzip',sep='\t')
        a_h.append(data)
        n=x.split(".")[0]
        d_h[n]=len(a_h)-1
    return a_h,d_h
