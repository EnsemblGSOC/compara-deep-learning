import sys
import os
from req_data import get_data_file,download_data
from read_data import read_data_genome,read_data_homology
from process_data import list_dict_genomes

def get_data_genome(arg,dir,a,d,ld,ldg):
    if arg[0]=='-d':
        if arg[4]=="-r":
            c=0
        else:
            return a,d,ld,ldg
    elif arg[0]=='-f':
        get_data_file(arg[1],dir)
    elif arg[0]=="-nd":
        return ld,ldg,a,d

    if arg[4]=="-r":
        a,d=read_data_genome(dir,a,d)
        assert(len(a)==len(d))

        ld,ldg=list_dict_genomes(a,d)

        assert(len(ld)==len(ldg))

        for i in range(len(ld)):
            assert(len(ld[i])==len(ldg[i]))

    return ld,ldg,a,d

def get_data_homology(arg,dir,a_h,d_h):
    if arg[2]=="-l":
        if not os.path.exists(dir):
            os.mkdir(dir)
        download_data(arg[3],dir)
    elif arg[2]=="-f":
        get_data_file(arg[3],dir)
    elif arg[2]=="-d":
        if arg[4]=="-r":
            c=0
        else:
            return a_h,d_h
    elif arg[2]=="-nd":
        return a_h,d_h

    if arg[4]=="-r":
        a_h,d_h=read_data_homology(dir,a_h,d_h)
        assert(len(a_h)==len(d_h))

    return a_h,d_h
