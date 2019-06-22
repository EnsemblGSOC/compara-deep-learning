import pandas as pd 
import requests 
import sys
import pickle
from get_data import get_data_genome

arg=sys.argv
arg=arg[1:]

if len(arg)!=3:
    print("No. of arguments more or less. Please check")
    sys.exit(1)

dir_g="data"
cmap,cimap,ld,ldg,a,d=get_data_genome(arg,dir_g)

data=dict(cmap=cmap,cimap=cimap,ld=ld,ldg=ldg,a=a,d=d)

with open("genome_maps","wb") as file:
    pickle.dump(data,file)

print("Genome Maps Created Successfully.")