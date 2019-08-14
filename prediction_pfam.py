import json
import gc
import pandas as pd 
import numpy as np 
import pickle
import tensorflow as tf
import sys
import os
import progressbar
from neighbor_genes import read_genome_maps
from process_data import create_data_homology_ls
from read_get_gene_seq import read_gene_sequences
from access_data_rest import update_rest,update_rest_protein
from threads import Procerssrunner
from prepare_synteny_matrix import read_data_synteny
from tree_data import create_tree_data
from process_data import create_map_list
from pfam_parser import pfam_parse
from pfam_matrix import create_pfam_map,create_pfam_matrix

def read_database(fname):
    df=pd.read_csv(fname,sep="\t",header=None)
    label_dict=dict(ortholog_one2one=1,
                    other_paralog=0,
                    non_homolog=2,
                    ortholog_one2many=1,
                    ortholog_many2many=1,
                    within_species_paralog=0,
                    gene_split=4)
    label=[]
    for _,row in df.iterrows():
        label.append(label_dict[row[7]])
    df=df.assign(label=label)
    df=df.drop(7,axis=1)
    df=df.drop(0,axis=1)
    df.columns=["gene_stable_id","species","homology_gene_stable_id","homology_species","goc","wga","label"]
    return df

def select_data_by_length(df,st,end):
    try:
        if end<len(df):
            if st<end:
                df=df.loc[df.index.values[st:end]]
            else:
                raise ValueError()
    except:
        print("Making Predictions for the complete dataframe:)")
    print(len(df))
    return df

def create_synteny_features(a_h,d_h,n,a,d,ld,ldg,cmap,cimap,name):    
    lsy=create_data_homology_ls(a_h,d_h,n,a,d,ld,ldg,cmap,cimap,0)
    gene_sequences=read_gene_sequences(a_h,lsy,"geneseq","prediction_"+name)
    gene_sequences=update_rest(gene_sequences,"prediction_"+name)
    print("Gene Sequences Loaded.")
    return lsy,gene_sequences

def threadmaker(nop,df,lsy,gene_sequences,n,name):
    part=len(df)//nop
    pr=Procerssrunner()
    pr.start_processes(nop,df,gene_sequences,lsy,part,n,name)
    smg,sml,indexes=read_data_synteny(nop,name)
    sml=np.array(sml)
    smg=np.array(smg)
    indexes=np.array(indexes)
    return sml,smg,indexes

def get_prediction(smg,sml,pfam_matrices,indexes,bls,blhs,dis,dps,dphs,model_name,no_of_model,w):
    preds=np.zeros((len(smg),no_of_model))
    pfam_matrices=pfam_matrices.reshape((len(smg),7,7,1))
    for i in range(1,no_of_model+1):
        try:
            model=tf.train.import_meta_graph(model_name+'_v'+str(i)+'/model.ckpt.meta')
        except:
            print("Something wrong with the model.")
            continue
        with tf.Session() as sess:
            try:
                model.restore(sess,model_name+'_v'+str(i)+"/model.ckpt")
                graph = tf.get_default_graph()
                synmgt,synmlt,pfamt,blst,blhst,dpst,dphst,dist,lrt,yt=graph.get_collection("input_nodes")
                predictions=graph.get_tensor_by_name("Predictions/BiasAdd:0")
                print("Model Loaded Successfully :)")
            except:
                print(":(")
                sys.exit()

            fd={synmgt:smg,
                synmlt:sml,
                pfamt:pfam_matrices,
                blst:bls,
                blhst:blhs,
                dpst:dps.reshape((len(blhs),1)),
                dist:dis.reshape((len(blhs),1)),
                dphst:dphs.reshape((len(blhs),1))}
            preds_t_1=sess.run([predictions],feed_dict=fd)
            preds_t_1=np.array(preds_t_1)[0]
            fd={synmgt:smg.transpose((0,2,1,3)),
                synmlt:sml.transpose((0,2,1,3)),
                pfamt:pfam_matrices.transpose((0,2,1,3)),
                blst:blhs,
                blhst:bls,
                dpst:dphs.reshape((len(blhs),1)),
                dist:dis.reshape((len(blhs),1)),
                dphst:dps.reshape((len(blhs),1))}
            preds_t_2=sess.run([predictions],feed_dict=fd)
            preds_t_2=np.array(preds_t_2)[0]
            preds=preds+w[i-1]*(preds_t_1+preds_t_2)/2
        tf.reset_default_graph()    
    preds=np.argmax(preds,axis=1)
    print(preds.shape)
    return preds

def write_preds(fname,model_name,name,preds,index_dict,df):
    print("Writing predcitions to:","prediction_"+fname+"_"+model_name+"_"+name+"_multiple_pfam.txt")
    with open("prediction_"+fname+"_"+model_name+"_"+name+"_multiple_pfam.txt","w") as file:
        for index,row in progressbar.progressbar(df.iterrows()):
            file.write(row["gene_stable_id"])
            file.write("\t")
            file.write(row["homology_gene_stable_id"])
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

def main():
    arg=sys.argv
    arg=arg[1:]
    fname=arg[0]
    model_name=arg[1]
    no_of_model=int(arg[2])
    nop=int(arg[3])
    st=int(arg[4])
    end=int(arg[5])
    name=arg[6]
    pfam_fname=arg[7]
    weight=arg[8]
    if weight=="e":
        w=[1]*no_of_model
    else:
        w=[]
        for i in range(no_of_model):
            w.append(float(arg[9+i]))
    df=read_database(fname)
    df=select_data_by_length(df,st,end)
    n=3

    if os.path.exists("prediction_data/data_"+fname):

        with open("prediction_data/data_"+fname,"rb") as file:
            save_dict=pickle.load(file)
        smg=save_dict["smg"]
        sml=save_dict["sml"]
        pfam_matrices=save_dict["pfam"]
        indexes=save_dict["indexes"]
        bls=save_dict["bls"]
        blhs=save_dict["blhs"]
        dis=save_dict["dis"]
        dps=save_dict["dps"]
        dphs=save_dict["dphs"]

    else:
        a,d,ld,ldg,cmap,cimap=read_genome_maps()#read the genome maps
        print("Genome Maps Loaded.")
        a_h=[df]
        d_h=["prediction"]

        lsy,gene_sequences=create_synteny_features(a_h,d_h,n,a,d,ld,ldg,cmap,cimap,name)
        sml,smg,indexes=threadmaker(nop,df,lsy,gene_sequences,n,name)
        a=""
        d=""
        ld=""
        ldg=""
        cmap=""
        cimap=""
        gene_sequences=""
        gc.collect()
        df_temp=df.loc[indexes]
        pfam_db=pd.read_hdf(pfam_fname.split(".")[0]+"_pfam_db.h5")
        with open(pfam_fname.split(".")[0]+"_pfam_map","rb") as file:
            pfam_map=pickle.load(file)
        pfam_matrices,indexes_pfam=create_pfam_matrix(df_temp,lsy,pfam_db,pfam_map)
        pfam_db=""
        gc.collect()
        bls,blhs,dis,dps,dphs=create_tree_data("species_tree.tree",df_temp)
        save_dict=dict(smg=smg,sml=sml,pfam=pfam_matrices,indexes=indexes,bls=bls,blhs=blhs,dis=dis,dps=dps,dphs=dphs)

        if not os.path.isdir("prediction_data"):
            os.mkdir("prediction_data")

        with open("prediction_data/data_"+fname,"wb") as file:
            pickle.dump(save_dict,file)
    
    index_dict=create_map_list(indexes)
    preds=get_prediction(smg,sml,pfam_matrices,indexes,bls,blhs,dis,dps,dphs,model_name,no_of_model,w)

    write_preds(fname,model_name,name,preds,index_dict,df)

if __name__=="__main__":
    main()