import pandas as pd
from threading import Thread
from multiprocessing import Process,Lock,Manager
import numpy as np
import edlib as ed
import pandas as pd
import progressbar
import time
from skbio.alignment import local_pairwise_align_ssw
from skbio import DNA,TabularMSA,RNA
import copy
import time
from save_data import write_data_synteny

class Thread_objects():
    def __init__(self,df_temp,gene_sequences,lsy,i,name):
        self.gene_sequences=copy.deepcopy(gene_sequences)
        self.lsy=copy.deepcopy(lsy)
        self.df=copy.deepcopy(df_temp)
        self.smg=[]
        self.sml=[]
        self.indexes=[]
        self.i=i
        self.start=0
        self.end=0
        self.name=name

    def create_synteny_matrix_mul(self,gene_seq,g1,g2,n):
        for gene in g1:
            if gene=="NULL_GENE":
                continue
            try:
                temp=gene_seq[gene]
            except:
                return np.zeros((n,n,2)),np.zeros((n,n,2))
        for gene in g2:
            if gene=="NULL_GENE":
                continue
            try:
                temp=gene_seq[gene]
            except:
                return np.zeros((n,n,2)),np.zeros((n,n,2))
        sm=np.zeros((n,n,2))
        sml=np.zeros((n,n,2))
        for i in range(n):
            if g1[i]=="NULL_GENE":
                continue
            if  gene_seq[g1[i]]=="":
                return np.zeros((n,n,2)),np.zeros((n,n,2))
            for j in range(n):
                if g2[j]=="NULL_GENE":
                    continue
                if  gene_seq[g2[j]]=="":
                    return np.zeros((n,n,2)),np.zeros((n,n,2))
                norm_len=max(len(gene_seq[g1[i]]),len(gene_seq[g2[j]]))
                try:
                    result = ed.align(gene_seq[g1[i]],gene_seq[g2[j]], mode="NW", task="distance")
                    sm[i][j][0]=result["editDistance"]/(norm_len)
                    result = ed.align(gene_seq[g1[i]],gene_seq[g2[j]][::-1], mode="NW", task="distance")
                    sm[i][j][1]=result["editDistance"]/(norm_len)
                    _,result,_=local_pairwise_align_ssw(DNA(gene_seq[g1[i]]),DNA(gene_seq[g2[j]]))
                    sml[i][j][0]=result/(norm_len)
                    _,result,_=local_pairwise_align_ssw(DNA(gene_seq[g1[i]]),DNA(gene_seq[g2[j]][::-1]))
                    sml[i][j][1]=result/(norm_len)
                except:
                    return np.zeros((n,n,2)),np.zeros((n,n,2))
        return sm,sml

    def synteny_matrix(self,gene_seq,hdf,lsy,n):
        t=0
        self.start=time.time()       
        for index,row in progressbar.progressbar(hdf.iterrows()):
            g1=str(row["gene_stable_id"])
            g2=str(row["homology_gene_stable_id"])
            x=[]
            y=[]
            t+=1
            try:
                temp=lsy[g1]
            except:
                continue
            try:
                temp=lsy[g2]
            except:
                continue
            for i in range(len(lsy[g1]['b'])-1,-1,-1):
                x.append(lsy[g1]['b'][i])
            x.append(g1)
            for k in lsy[g1]['f']:
                x.append(k)

            for i in range(len(lsy[g2]['b'])-1,-1,-1):
                y.append(lsy[g2]['b'][i])
            y.append(g2)
            for k in lsy[g2]['f']:
                y.append(k)
                
            assert(len(x)==len(y))
            assert(len(x)==(2*n+1))
            smgtemp,smltemp=self.create_synteny_matrix_mul(gene_seq,x,y,2*n+1)
            if np.all(smgtemp==0):
                continue  
            self.smg.append(smgtemp)
            self.sml.append(smltemp)
            self.indexes.append(index)
        self.end=time.time()
        print("Thread {} finished in {}s.".format(self.i+1,self.end-self.start))
        write_data_synteny(self.smg,self.sml,self.indexes,self.i,self.name)

class Procerssrunner():
    def __init__(self):
        self.thread_alive=[]
        self.obj_list=[]    

    def start_thread(self,obj,i,thread_alive,n,name):
        t=Process(target=obj.synteny_matrix,args=(obj.gene_sequences,obj.df,obj.lsy,n),name="Thread_"+str(i+1))
        print("Thread ",(i+1)," started for ",name,".")
        thread_alive.append(t)

    def start_processes(self,nop,df,gene_sequences,lsy,part,n,name):
        for i in range(nop):
            df_temp=df.loc[df.index.values[i*part:(i+1)*part]]
            obj=Thread_objects(df_temp,gene_sequences,lsy,i,name)
            print("Object Created")
            self.start_thread(obj,i,self.thread_alive,n,name)
            self.obj_list.append(obj)
        st=time.time()
        for t in self.thread_alive:
            t.start()
        #for t in self.thread_alive:
            #t.join()
        while(len(self.thread_alive)!=0):
            time.sleep(60)
            self.thread_alive=[t for t in self.thread_alive if t.is_alive()]
        end=time.time()
        print("Ending Processes")
        print("Time taken:{}s".format(end-st))

