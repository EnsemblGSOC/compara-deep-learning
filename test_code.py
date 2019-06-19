import sys
import numpy as np
import pandas as pd
import os
import gc
from get_data import get_data_homology,get_data_genome

arg=sys.argv
arg=arg[1:]

if len(arg)!=5:
    print("No. of arguments more or less. Please check")
    sys.exit(1)

dir_g="data"
cmap,cimap,ld,ldg,a,d=get_data_genome(arg,dir_g)

print(d)
print (cmap[0]['ENSG00000219481'])
print(cimap[0]['1'])