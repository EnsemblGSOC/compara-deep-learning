from ete3 import Tree
import numpy as np
import progressbar

def create_branch_length_padding(bl):
    maxlen=0
    for x in bl:
        if len(x)>maxlen:
            maxlen=len(x)

    for x in bl:
        for i in range(len(x),maxlen):
            x.append(0)

def create_tree_data(treename,df):
    t=Tree(treename)
    branch_lengths_s=[]
    branch_lengths_hs=[]
    dist=[]
    ns=[]
    nhs=[]
    for index,row in progressbar.progressbar(df.iterrows()):
        d=0
        x=row["species"]
        y=row["homology_species"]
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
