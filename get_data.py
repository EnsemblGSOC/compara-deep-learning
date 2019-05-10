import sys
import os
from req_data import get_data_file
from read_data import read_data_genome
from process_data import list_dict_genomes

arg=sys.argv
arg=arg[1:]

a=[]
d={}

feature_name='gene'

if arg[0]=='-d':
    a,d=read_data_genome(arg[1],a,d)
elif arg[0]=='-f':
    dir_name=get_data_file(arg[1])
    a,d=read_data_genome(dir_name,a,d)

assert(len(a)==len(d))

ld,ldg=list_dict_genomes(a,d)

assert(len(ld)==len(ldg))

for i in range(len(ld)):
    assert(len(ld[i])==len(ldg[i]))
