import sys
import numpy as np
import pandas as pd
import json
import os
import gc
from get_data import get_data_homology,get_data_genome
from process_data import create_data_homology_ls

if not os.path.exists("processed"):
    os.mkdir("processed")

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
cmap,cimap,ld,ldg,a,d=get_data_genome(arg,dir_g)

dir_hom="data_homology"
a_h,d_h=get_data_homology(arg,dir_hom)

if arg[-1]=="-d":
    sys.exit(1)

print("Data Read")

n=3 #no. of numbers neighbors
save_after=5 #to save data after n steps

lsy=create_data_homology_ls(a_h,d_h,n,a,d,ld,ldg,cmap,cimap,save_after,enable_break,1)
print(len(lsy))

print("Neighbor Genes Updated Successfully")



