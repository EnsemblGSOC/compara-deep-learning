import tensorflow as tf 
import numpy as np 
import pandas as pd 
import requests 
import sys
import pickle
from get_data import get_data_genome
from process_data import get_nearest_neighbors
from create_synteny_matrix_v2 import create_synteny_matrix_mul
from ete3 import Tree
import traceback

"""
arg=sys.argv
arg=arg[1:]

enable_break=0

if arg[-1]=="-test":
    enable_break=1

arg=arg[:-1]

if len(arg)!=5:
    print("No. of arguments more or less. Please check")
    sys.exit(1)

dir_g="data"
cmap,cimap,ld,ldg,a,d=get_data_genome(arg,dir_g)"""

data={}
with open("genome_maps","rb") as file:
    data=pickle.load(file)
cmap=data["cmap"]
cimap=data["cimap"]
ld=data["ld"]
ldg=data["ldg"]
a=data["a"]
d=data["d"]

n=3

server = "https://rest.ensembl.org"
def get_gene_data(gid):
    empty={}
    ext = "/lookup/id/"+gid+"?" 
    try:
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"}) 
        if not r.ok:
            r.raise_for_status()
    except:
        return empty
    decoded = r.json()
    return dict(decoded)

def create_branch_length_padding(bl):
    maxlen=29
    for x in bl:
        for i in range(len(x),maxlen):
            x.append(0)

def create_tree_data(x,y):
    t=Tree("species_tree.tree")
    branch_lengths_s=[]
    branch_lengths_hs=[]
    dist=[]
    ns=[]
    nhs=[]
    d=0
    bl=[]
    c=0
    mca=t.get_common_ancestor(x,y)
    node=t&x
    while node.up!=mca:
        d+=node.dist
        bl.append(node.dist)
        node=node.up
        c+=1
    ns.append(c)
    c=0
    branch_lengths_s.append(bl)
    bl=[]
    node=t&y
    while node.up!=mca:
        d+=node.dist
        bl.append(node.dist)
        node=node.up
        c+=1
    nhs.append(c)
    branch_lengths_hs.append(bl)
    dist.append(d)
    create_branch_length_padding(branch_lengths_s)
    create_branch_length_padding(branch_lengths_hs)
    return np.array(branch_lengths_s),np.array(branch_lengths_hs),np.array(dist),np.array(ns),np.array(nhs)

def filter_gene_details(gid,n,a,d,ld,ldg,cmap,cimap):
    gd=get_gene_data(gid)
    if len(gd)==0 or gd["biotype"]!="protein_coding":
        print("Gene Id is Incorrect!!!. Try again.")
        return list()
    sp=gd["species"]
    print(sp)
    try:
        t2=d[sp.capitalize()]
    except:
        print("Gene Species does not exist in gene maps. Try another species.")
        return list(),sp
    nl,nr=get_nearest_neighbors(gid,sp,n,a,d,ld,ldg,cmap,cimap)
    if len(nl)==0 or len(nr)==0:
        print("Gene is not present in gtf file. Please Try another gene.")
        return list(),sp
    print("Found nearest neighbors")
    x=[]
    for i in range(len(nl)-1,-1,-1):
        x.append(nl[i])
    x.append(gid)
    for k in nr:
        x.append(k)
    assert(len(x)==2*n+1)
    return x,sp

def get_label(l):
    if l==0:
        return "Paralogous"
    if l==1:
        return "Orthologous"
    return "Non Homologous"

def softmax(arr):
    arr=arr-np.amax(arr)
    arr=np.exp(arr)
    arr=arr/np.sum(arr)
    return arr

geneseq={}
try:
    model=tf.train.import_meta_graph('saved_models/model.ckpt.meta')
except:
    print("Something wrong with the model.")
    sys.exit(1)
with tf.Session() as sess:
    try:
        model.restore(sess,"saved_models/model.ckpt")
        graph = tf.get_default_graph()
        synmgt,synmlt,blst,blhst,dpst,dphst,dist,lrt,yt=graph.get_collection("input_nodes")
        predictions=graph.get_tensor_by_name("Predictions/BiasAdd:0")
        print("Model Loaded Successfully :)")
    except:
        print(":(")
        sys.exit()
    while(1):
        try:
            ch=input("Do you want to enter a gene id [y/n]:")
            if ch=='n':
                break
            g1=input("Enter the first gene id:")
            g1,sp1=filter_gene_details(g1,3,a,d,ld,ldg,cmap,cimap)
            if len(g1)==0:
                continue
            g2=input("Enter the second gene id:")
            g2,sp2=filter_gene_details(g2,3,a,d,ld,ldg,cmap,cimap)
            if len(g2)==0:
                continue
            smg,sml=create_synteny_matrix_mul(geneseq,g1,g2,2*n+1)
            bls,blhs,dis,dps,dphs=create_tree_data(sp1,sp2)
            fd={synmgt:smg.reshape((1,7,7,2)),
                synmlt:sml.reshape((1,7,7,2)),
                blst:bls.reshape((1,29)),
                blhst:blhs.reshape((1,29)),
                dpst:dps.reshape((1,1)),
                dist:dis.reshape((1,1)),
                dphst:dphs.reshape((1,1))}
            preds=sess.run([predictions],feed_dict=fd)
            label=np.argmax(preds)
            print(get_label(label))
            print(softmax(preds))
            
        except Exception as e:
            print("Some Error Was There Try Again:(",e)
            traceback.print_exc()
            continue