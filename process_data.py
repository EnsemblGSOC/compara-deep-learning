import pandas
import gc
import numpy as np

def get_nearest_neighbors(g,gs,n,a,d,ld,ldg):
    ne=[]
    nr=[]
    gi=d[gs.capitalize()]
    sldf=a[gi]
    sld=ld[gi]
    sldg=ldg[gi]
    if g not in sldg:
        return ne
    i=sldg[g]
    #get the -n neighbors
    start=int(sldf.iloc[i,[3]])
    for j in range(n):
        #select the column
        itemp=0
        end=sldf.iloc[:,4]
        end=np.array(end)
        assert(len(end)==len(sld))
        end=end-start
        end_s=np.argsort(end)
        for k in end_s:
            if end[k]<0 and end[k+1]>=0:
                itemp=k
                break
        ne.append(sld[itemp])
        start=int(sldf.iloc[itemp,[3]])
        #print(start)
    #get the +n neighbors
    end=int(sldf.iloc[i,[4]])
    for j in range(n):
        itemp=0
        start=sldf.iloc[:,3]
        start=np.array(start)
        start=start-end
        start_s=np.argsort(start)
        for k in start_s:
            if start[k]>0:
                itemp=k
                break
        nr.append(sld[itemp])
        end=int(sldf.iloc[itemp,[4]])

    return ne,nr

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

def create_data_homology_ls(a_h,d_h,n,a,d,ld,ldg):
    lsy={} #dictionary which stores +/- n genes of the given gene by id. Each key is a gene id which corresponds to the one in center.
    lcmap={} #dictionary which stores the gene pairs already considered
    t=0
    for df in a_h:
        for index,row in df.iterrows():
            x=row["gene_stable_id"]
            y=row["homology_gene_stable_id"]
            xs=row["species"]
            ys=row["homology_species"]
            if x+y in lcmap or y+x in lcmap:
                continue
            if x not in lsy:
                xl,xr=get_nearest_neighbors(x,xs,n,a,d,ld,ldg)
                lsy[x]=dict(b=xl,f=xr)
            if y not in lsy:
                yarr=[]
                yl,yr=get_nearest_neighbors(y,ys,n,a,d,ld,ldg)
                lsy[y]=dict(b=yl,f=yr)
            lcmap[x+y]=1


    return lsy,lcmap
