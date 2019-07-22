import json
import gc
import pandas as pd 
import numpy as np 
import pickle
import tensorflow as tf
import sys
import time
from test_prepare_functions import create_data_homology_ls,read_gene_sequences,create_tree_data,read_data,update_rest
from process_data import create_map_list
from threading import Thread,Lock
import traceback
from threads import Procerssrunner

def main():
    arg=sys.argv
    fname=arg[-6]
    model_name=arg[-5]
    n_of_t=int(arg[-4])
    st=int(arg[-3])
    end=int(arg[-2])
    name=arg[-1]

    df=pd.read_csv(fname,sep="\t",header=None)
    label_dict=dict(ortholog_one2one="Orthologs",
                    other_paralog="Paralogs",
                    ortholog_one2many="Orthologs",
                    ortholog_many2many="Orthologs",
                    within_species_paralog="Paralogs",
                    gene_split="Gene Split",
                    non_homolog="non_homolog")
    label_dict_2=dict(ortholog_one2one=1,other_paralog=0,non_homolog=2,ortholog_one2many=1,ortholog_many2many=1,within_species_paralog=0,gene_split=4)
    label_2=[]
    label_1=[]
    for _,row in df.iterrows():
        label_1.append(label_dict_2[row[7]])
        label_2.append(label_dict[row[7]])
    df=df.assign(label=label_1)
    df[7]=label_2
    try:
        if end<len(df):
            if st<end:
                df=df.loc[df.index.values[st:end]]
            else:
                raise ValueError()
    except:
        print("Making Predictions for the complete dataframe:)")
    print(len(df))
    
    data={}    
    with open("genome_maps","rb") as file:
        data=pickle.load(file)
    cmap=data["cmap"]
    cimap=data["cimap"]
    ld=data["ld"]
    ldg=data["ldg"]
    a=data["a"]
    d=data["d"]
    lsy=create_data_homology_ls(df,3,a,d,ld,ldg,cmap,cimap)
    print(len(lsy))
    data={}
    a=[]
    d=[]
    ldg=[]
    ld=[]
    cmap=[]
    cimap=[]
    gc.collect()
    
    gene_sequences=read_gene_sequences(df,lsy,"geneseq","gene_sequences",name)
    print("Gene Sequences Loaded")

    print("Going to update not found sequences:")
    gene_sequences=update_rest(gene_sequences,name)
    
    gc.collect()
    part=len(df)//n_of_t
    """
    lsy={}
    gene_sequences={}
    part=1"""
    return lsy,gene_sequences,df,model_name,fname,n_of_t,part,name


if __name__=='__main__':
    n=3
    lsy,gene_sequences,df,model_name,fname,n_of_t,part,name=main()
    pr=Procerssrunner()
    pr.start_processes(n_of_t,df,gene_sequences,lsy,part,n,name)
    smg,sml,indexes=read_data(n_of_t,name)
    sml=np.array(sml)
    smg=np.array(smg)
    indexes=np.array(indexes)
    print(indexes.shape)
    df_temp=df.loc[indexes]
    bls,blhs,dis,dps,dphs=create_tree_data("species_tree.tree",df_temp)

    assert(len(bls)==len(indexes))
    assert(len(blhs)==len(smg))
    assert(len(df_temp)==len(dps))

    print("Data Prepared.")
    index_dict=create_map_list(indexes)    
    preds=np.zeros((len(smg),3))
    w=[0.86,0.8,0.06,0.0]
    for i in range(1,4):
        try:
            model=tf.train.import_meta_graph(model_name+'_v'+str(i)+'/model.ckpt.meta')
        except:
            print("Something wrong with the model.")
            continue
        with tf.Session() as sess:
            try:
                model.restore(sess,model_name+'_v'+str(i)+"/model.ckpt")
                graph = tf.get_default_graph()
                synmgt,synmlt,blst,blhst,dpst,dphst,dist,lrt,yt=graph.get_collection("input_nodes")
                predictions=graph.get_tensor_by_name("Predictions/BiasAdd:0")
                print("Model Loaded Successfully :)")
            except:
                print(":(")
                sys.exit()

            fd={synmgt:smg,
                synmlt:sml,
                blst:bls,
                blhst:blhs,
                dpst:dps.reshape((len(blhs),1)),
                dist:dis.reshape((len(blhs),1)),
                dphst:dphs.reshape((len(blhs),1))}
            preds_t_1=sess.run([predictions],feed_dict=fd)
            preds_t_1=np.array(preds_t_1)[0]
            fd={synmgt:smg.transpose((0,2,1,3)),
                synmlt:sml.transpose((0,2,1,3)),
                blst:blhs,
                blhst:bls,
                dpst:dphs.reshape((len(blhs),1)),
                dist:dis.reshape((len(blhs),1)),
                dphst:dps.reshape((len(blhs),1))}
            preds_t_2=sess.run([predictions],feed_dict=fd)
            preds_t_2=np.array(preds_t_2)[0]
            preds=preds+w[i-1]*(preds_t_1+preds_t_2)/2
        tf.reset_default_graph()
    
    #preds=preds/6
    preds=np.argmax(preds,axis=1)
    print(preds.shape)
        
    with open("prediction_"+fname+"_"+model_name+"_"+name+"_multiple.txt","w") as file:
        for index,row in df.iterrows():
            file.write(str(row[0]))
            file.write("\t")
            file.write(row[1])
            file.write("\t")
            file.write(row[3])
            file.write("\t")
            file.write(str(row["label"]))
            file.write("\t")
            if index in index_dict:
                file.write(str(preds[index_dict[index]]))
                file.write("\t")
                if preds[index_dict[index]]==row["label"]:
                    file.write(str(1))
                else:
                    file.write(str(0))
            else:
                file.write("Error")
                file.write("\t")
                file.write("NaN")
            file.write("\n")