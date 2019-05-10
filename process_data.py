import pandas
import gc


ls=[]
ld=[]
def list_dict_genomes(a,n):
    for x in a:
        ldg={}
        uc=list(x["gene_id"])
        for i in range(len(uc)):
            ldg[uc[i]]=i
        ls.append(uc)
        ld.append(ldg)
    return ls,ld
