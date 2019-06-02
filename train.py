import tensorflow as tf
import numpy as np

from model import create_model

def train(train_synteny_matrices_global,train_synteny_matrices_local,train_branch_length_species,train_branch_length_homology_species,train_mean_gene_length,train_dist_p_s,train_dist_p_hs,train_distance,train_labels):
    graph,saver=create_model()
    synmg,synml,bls,blhs,gl,dps,dphs,dis,lr,y=graph.get_collection("input_nodes")
    loss,t_op,accuracy,init=graph.get_collection("output_nodes")
    with tf.Session(graph=graph) as sess:
        sess.run(init)
        batch_size=64
        num_epochs=30
        learn=0.001
        for j in range(num_epochs):
            for i in range(50000//batch_size):
                feed_dict={
                    synmg:train_synteny_matrices_global[i*batch_size:(i+1)*batch_size],
                    synml:train_synteny_matrices_local[i*batch_size:(i+1)*batch_size],
                    bls:train_branch_length_species[i*batch_size:(i+1)*batch_size],
                    blhs:train_branch_length_homology_species[i*batch_size:(i+1)*batch_size],
                    gl:train_mean_gene_length[i*batch_size:(i+1)*batch_size].reshape((batch_size,1)),
                    dps:train_dist_p_s[i*batch_size:(i+1)*batch_size].reshape((batch_size,1)),
                    dphs:train_dist_p_hs[i*batch_size:(i+1)*batch_size].reshape((batch_size,1)),
                    dis:train_distance[i*batch_size:(i+1)*batch_size].reshape((batch_size,1)),
                    lr:learn,
                    y:train_labels[i*batch_size:(i+1)*batch_size]
                }
                sess.run(t_op,feed_dict=feed_dict)

            test_dict={
                    synmg:train_synteny_matrices_global[-15000:],
                    synml:train_synteny_matrices_local[-15000:],
                    bls:train_branch_length_species[-15000:],
                    blhs:train_branch_length_homology_species[-15000:],
                    gl:train_mean_gene_length[-15000:].reshape((15000,1)),
                    dps:train_dist_p_s[-15000:].reshape((15000,1)),
                    dphs:train_dist_p_hs[-15000:].reshape((15000,1)),
                    dis:train_distance[-15000:].reshape((15000,1)),
                    y:train_labels[-15000:]
                    }
            accuracy_test,loss_test=sess.run([accuracy,loss],feed_dict=test_dict)
            print("Epoch:{} Test Accuracy:{} Test Loss:{}".format(j+1,accuracy_test*100,loss_test))
            learn*=0.97
        saver.save(sess,"saved_models/model.ckpt")
