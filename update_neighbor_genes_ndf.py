import sys
import numpy as np
import pandas as pd
import json
import os
import gc
import pickle
from get_data import get_data_homology,get_data_genome
from process_data import create_data_homology_ls

data={}
with open("genome_maps","rb") as file:
    data=pickle.load(file)
cmap=data["cmap"]
cimap=data["cimap"]
ld=data["ld"]
ldg=data["ldg"]
a=data["a"]
d=data["d"]

df=pd.read_hdf("negative_dataset.h5",key="ndf")
a_h=[]
d_h=[]
a_h.append(df)
d_h.append("negative_dataset")
print("Data Read")

n=3 #no. of numbers neighbors
save_after=5 #to save data after n steps

lsy=create_data_homology_ls(a_h,d_h,n,a,d,ld,ldg,cmap,cimap,save_after,0)
print(len(lsy))
print("Neighbor Genes Updated Successfully")