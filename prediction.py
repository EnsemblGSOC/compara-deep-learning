import json
import gc
import pandas as pd 
import numpy as np 
import pickle
import tensorflow as tf
import sys
from test_prepare_functions import create_data_homology_ls,read_gene_sequences,synteny_matrix,create_tree_data
from access_data_rest import update_rest
from process_data import create_map_list

df=pd.read_csv("input_real_predictions_random_20K.txt",sep="\t",header=None)
label_dict=dict(ortholog_one2one="Orthologs",
                other_paralog="Paralogs",
                ortholog_one2many="Orthologs",
                ortholog_many2many="Orthologs",
                within_species_paralog="Paralogs",
                gene_split="Gene Split")
label_dict_2=dict(ortholog_one2one=1,other_paralog=0,ortholog_one2many=1,ortholog_many2many=1,within_species_paralog=0,gene_split=4)
label_2=[]
label_1=[]
for index,row in df.iterrows():
    label_1.append(label_dict_2[row[7]])
    label_2.append(label_dict[row[7]])

df=df.assign(label=label_1)
df[7]=label_2
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
with open("ng","wb")as file:
    pickle.dump(lsy,file)

gene_sequences=read_gene_sequences(df,lsy,"geneseq","gene_sequences")
print("Gene Sequences Loaded")

print("Going to update not found sequences:")
gene_sequences=update_rest(gene_sequences)
with open("gs","wb")as file:
    pickle.dump(gene_sequences,file)

"""with open("ng","rb") as file:
    lsy=pickle.load(file)
with open("gs","rb") as file:
    gene_sequences=pickle.load(file)"""

n=3
smg,sml,indexes=synteny_matrix(gene_sequences,df,lsy,n,0,list(),list(),list())
print(len(indexes))

df_temp=df.loc[indexes]
bls,blhs,dis,dps,dphs=create_tree_data("species_tree.tree",df_temp)

assert(len(bls)==len(indexes))
assert(len(blhs)==len(smg))
assert(len(df_temp)==len(dps))

print("Data Prepared.")

index_dict=create_map_list(indexes)

try:
    model=tf.train.import_meta_graph('saved_models/model.ckpt.meta')
except:
    print("Something wrong with the model.")
    sys.exit(1)

with tf.Session() as sess:
    try:
        model.restore(sess,"saved_models/model.ckpt")
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
    preds=sess.run([predictions],feed_dict=fd)
    preds=np.array(preds)[0]
    preds=np.argmax(preds,axis=1)
    print(preds.shape)
    
with open("prediction.txt","w") as file:
    for index,row in df.iterrows():
        file.write(str(row[0]))
        file.write("\t")
        file.write(row[1])
        file.write("\t")
        file.write(row[2])
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





