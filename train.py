import pickle
import numpy as np
import sys
import tensorflow as tf
from model import create_model

def create_branch_length_padding(bl):
    maxlen=0
    for x in bl:
        if len(x)>maxlen:
            maxlen=len(x)

    for x in range(len(bl)):
        temp=bl[x]
        for i in range(len(temp),maxlen):
            temp=np.append(temp,[0])
        bl[x]=temp
    return bl

def train(train_synteny_matrices_global,train_synteny_matrices_local,train_pfam_matrices,train_branch_length_species,train_branch_length_homology_species,train_dist_p_s,train_dist_p_hs,train_distance,train_labels,v,num_epochs,learning_rate,decay,size_train,batch_size,model_name):
    print("Going to train model {} for:\n Batch Size:{} \n Learning Rate:{} \n Decay:{}\n On {} Samples".format(v,batch_size,learning_rate,decay,len(train_synteny_matrices_global)))
    graph,saver=create_model()
    synmg,synml,pfam,bls,blhs,dps,dphs,dis,lr,y=graph.get_collection("input_nodes")
    loss,t_op,accuracy,init,summary=graph.get_collection("output_nodes")
    with tf.Session(graph=graph) as sess:
        writer = tf.summary.FileWriter('./'+model_name+'_v'+str(v), sess.graph)
        sess.run(init)
        learn=learning_rate
        for j in range(num_epochs):
            for i in range(size_train//batch_size):
                feed_dict={
                    synmg:train_synteny_matrices_global[i*batch_size:(i+1)*batch_size],
                    synml:train_synteny_matrices_local[i*batch_size:(i+1)*batch_size],
                    pfam:train_pfam_matrices[i*batch_size:(i+1)*batch_size].reshape((batch_size,7,7,1)),
                    bls:train_branch_length_species[i*batch_size:(i+1)*batch_size],
                    blhs:train_branch_length_homology_species[i*batch_size:(i+1)*batch_size],
                    dps:train_dist_p_s[i*batch_size:(i+1)*batch_size].reshape((batch_size,1)),
                    dphs:train_dist_p_hs[i*batch_size:(i+1)*batch_size].reshape((batch_size,1)),
                    dis:train_distance[i*batch_size:(i+1)*batch_size].reshape((batch_size,1)),
                    lr:learn,
                    y:train_labels[i*batch_size:(i+1)*batch_size]
                }
                sess.run(t_op,feed_dict=feed_dict)

            test_dict={
                    synmg:train_synteny_matrices_global[size_train:],
                    synml:train_synteny_matrices_local[size_train:],
                    pfam:train_pfam_matrices[size_train:].reshape((len(train_synteny_matrices_global)-size_train,7,7,1)),
                    bls:train_branch_length_species[size_train:],
                    blhs:train_branch_length_homology_species[size_train:],
                    dps:train_dist_p_s[size_train:].reshape((len(train_synteny_matrices_global)-size_train,1)),
                    dphs:train_dist_p_hs[size_train:].reshape((len(train_synteny_matrices_global)-size_train,1)),
                    dis:train_distance[size_train:].reshape((len(train_synteny_matrices_global)-size_train,1)),
                    y:train_labels[size_train:],
                    lr:learn
                    }
            accuracy_test,loss_test,summary_write=sess.run([accuracy,loss,summary],feed_dict=test_dict)
            writer.add_summary(summary_write,i+1)
            print("Epoch:{} Test Accuracy:{} Test Loss:{}".format(j+1,accuracy_test*100,loss_test))
            learn*=decay
        saver.save(sess,model_name+"_v"+str(v)+"/model.ckpt")
        writer.close()

def read_positive(len_p,bls,blhs,dis,dps,dphs,sml,smg,pfam,label):
    rowsh=[]
    with open("dataset","rb") as file:
        rowsh=pickle.load(file)
    shi=np.random.permutation(len(rowsh))
    rows_shuffled=[]
    for i in range(len(shi)):
        rows_shuffled.append(rowsh[shi[i]])  
    rowsh=rows_shuffled
    spco={}
    spcp={}
    for row in rowsh:
        if row["species"] not in spco:
            spco[row["species"]]=0
            spcp[row["species"]]=0
           
    maxcount_o=int((len_p)*0.3/14)
    maxcount_p=int((len_p)*0.7/14)
    for row in rowsh:
        if row["label"]==2:
            continue
        if row["label"]==1 and spco[row["species"]]>maxcount_o:
            continue
        if row["label"]==0 and spcp[row["species"]]>maxcount_p:
            continue
        bls.append(np.array(row["bls"]))
        blhs.append(np.array(row["blhs"]))
        dis.append(row["dis"])
        dps.append(row["dps"])
        dphs.append(row["dphs"])
        sml.append(row["local_alignment_matrix"])
        smg.append(row["global_alignment_matrix"])
        pfam.append(row["pfam_matrix"])
        label.append(row["label"])
        if row["label"]==1:
            spco[row["species"]]+=1
        if row["label"]==0:
            spcp[row["species"]]+=1
            
            
def read_negative(len_n,bls,blhs,dis,dps,dphs,sml,smg,pfam,label):
    rows=[]
    with open("dataset","rb") as file:
            rows=pickle.load(file)
    rows=[row for row in rows if row["label"]==2]
    shi=np.random.permutation(len(rows))
    rows_shuffled=[]
    for i in range(len(shi)):
        rows_shuffled.append(rows[shi[i]])
    rows=rows_shuffled
    rows=rows[:len_n]
    for row in rows:
        bls.append(np.array(row["bls"]))
        blhs.append(np.array(row["blhs"]))
        dis.append(row["dis"])
        dps.append(row["dps"])
        dphs.append(row["dphs"])
        sml.append(row["local_alignment_matrix"])
        smg.append(row["global_alignment_matrix"])
        label.append(row["label"])
        pfam.append(row["pfam_matrix"])

def train_models(model_name,start,end,num_epochs,learn_rate,decay,size_train,batch_size):
    k=1
    for i in range(start//10,end//10+1):
        bls=[]
        blhs=[]
        dis=[]
        dps=[]
        dphs=[]
        sml=[]
        smg=[]
        pfam=[]
        label=[]
        portion=float(i/10)
        len_n=int(size_train*portion)
        len_p=int(size_train*(1-portion))
        read_positive(len_p,bls,blhs,dis,dps,dphs,sml,smg,pfam,label)
        read_negative(len_n,bls,blhs,dis,dps,dphs,sml,smg,pfam,label)
        bls=create_branch_length_padding(bls)
        blhs=create_branch_length_padding(blhs)
        bls=np.array(bls)
        print(bls.shape)
        blhs=np.array(blhs)
        print(blhs.shape)
        dis=np.array(dis)
        print(dis.shape)
        dps=np.array(dps)
        print(dps.shape)
        dphs=np.array(dphs)
        print(dphs.shape)
        sml=np.array(sml)
        print(sml.shape)
        smg=np.array(smg)
        print(smg.shape)
        pfam=np.array(pfam)
        print(pfam.shape)
        label=np.array(label)
        print(label.shape)
        shi=np.random.permutation(len(label))
        labels=label[shi]
        bls=bls[shi]
        blhs=blhs[shi]
        dis=dis[shi]
        dps=dps[shi]
        dphs=dphs[shi]
        sml=sml[shi]
        smg=smg[shi]
        pfam=pfam[shi]
        train(smg,sml,pfam,bls,blhs,dps,dphs,dis,labels,k,num_epochs,learn_rate,decay,int(0.9*size_train),batch_size,model_name)
        k+=1

def main():
    arg=sys.argv
    model_name=arg[-8]
    start_p=int(arg[-7])
    end_p=int(arg[-6])
    num_epochs=int(arg[-5])
    learn_rate=float(arg[-4])
    decay=float(arg[-3])
    size_train=float(arg[-2])
    batch_size=int(arg[-1])
    train_models(model_name,start_p,end_p,num_epochs,learn_rate,decay,size_train,batch_size)
        

if __name__=="__main__":
    main()


