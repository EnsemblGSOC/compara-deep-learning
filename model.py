dim=10
n=3
fl=2*n+1
import tensorflow as tf

def create_model():
    tf.reset_default_graph()
    g=tf.Graph()
    with g.as_default():
        synmg=tf.placeholder(dtype=tf.float64,shape=(None,2*n+1,2*n+1,2),name="Synteny_matrix_placeholder_Global")
        synml=tf.placeholder(dtype=tf.float64,shape=(None,2*n+1,2*n+1,2),name="Synteny_matrix_placeholder_Local")
        bls=tf.placeholder(dtype=tf.float64,shape=(None,21),name="Species_Branch_Length_Placeholder")
        blhs=tf.placeholder(dtype=tf.float64,shape=(None,21),name="Homology_Species_Branch_Length_Placeholder")
        gl=tf.placeholder(dtype=tf.float64,shape=(None,1),name="Mean_gene_length")
        dps=tf.placeholder(dtype=tf.float64,shape=(None,1),name="mca_species_distance")
        dphs=tf.placeholder(dtype=tf.float64,shape=(None,1),name="mca_homology_species_distance")
        dis=tf.placeholder(dtype=tf.float64,shape=(None,1),name="total_distance")
        lr=tf.placeholder_with_default(0.1,(),"learning_rate")
        y=tf.placeholder(dtype=tf.int32,shape=(None),name="labels")

        #reshape the variables
        #gl=tf.expand_dims(gl,-1)
        #dps=tf.expand_dims(dps,-1)
        #dphs=tf.expand_dims(dphs,-1)
        #dis=tf.expand_dims(dis,-1)

        x=tf.concat([gl,dps-dphs,dis],1,name="Create_train_vector")

        print(synmg,"\n",synml,"\n",bls,"\n",blhs,"\n",gl,"\n",dps,"\n",dphs,"\n",dis,"\n",x)

        reg_l2=tf.contrib.layers.l2_regularizer(0.001)
        reg_l1 = tf.contrib.layers.l1_regularizer(scale=0.005, scope=None)

        def get_variable_by_shape(shape,name):
            f=tf.get_variable(name,shape=shape,initializer=tf.glorot_uniform_initializer(),dtype=tf.float64,regularizer=reg_l2)
            return f

        def create_synteny_aligner(name,synm):
            with tf.variable_scope(name+"Synteny_Aligner",reuse=tf.AUTO_REUSE):

                fconv=get_variable_by_shape((2,2,2,dim),"fconv")
                conv=tf.nn.conv2d(synm,fconv,(1,1,1,1),padding="VALID",name="Conv_aligner")

                fconv_1=get_variable_by_shape((2,2,dim,dim*2),"fconv_1")
                conv_1=tf.nn.conv2d(conv,fconv_1,(1,1,1,1),padding="VALID",name="Conv_aligner_1")
                print(conv)
                print(conv_1)
                
                
                #fxconv=get_variable_by_shape((1,5,2,dim),"fxconv")
                fxconv=get_variable_by_shape((fl,2,dim*2),"fxconv")
                x_conv=tf.reshape(synm,(-1,fl*fl,2))
                x_conv=tf.nn.conv1d(x_conv,fxconv,stride=fl,padding="SAME",name="row_aligner")
                print(x_conv)
                
                
                #fyconv=get_variable_by_shape((5,1,2,dim),"fyconv")
                #y_conv=tf.nn.conv2d(synm,fyconv,(1,1,1,1),padding="SAME",name="column_aligner")
                y_conv=tf.reshape(tf.transpose(synm,(0,2,1,3)),(-1,fl*fl,2))
                print(y_conv)
                fyconv=get_variable_by_shape((fl,2,dim*2),"fyconv")
                y_conv=tf.nn.conv1d(y_conv,fyconv,stride=fl,padding="SAME",name="column_aligner")
                print(y_conv)
                
                
                wconv=get_variable_by_shape((fl,fl,2,dim*2*10),"wconv")
                w_conv=tf.nn.conv2d(synm,wconv,(1,1,1,1),padding="VALID",name="Global_Aligner_1")
                print(w_conv)
                
                conv_1=tf.reshape(conv_1,(-1,25,dim*2))
                print(conv)
                x_conv=tf.reshape(x_conv,(-1,fl,dim*2))
                print(x_conv)
                y_conv=tf.reshape(y_conv,(-1,fl,dim*2))
                print(y_conv)
                w_conv=tf.reshape(w_conv,(-1,10,dim*2))
                print(wconv)
                conv_final=tf.concat([conv_1,x_conv,y_conv,w_conv],1,name="Concatenate_All_Alignments")
                print(conv_final)
            return conv_final

        conv_final_g=create_synteny_aligner("Global",synmg)
        conv_final_l=create_synteny_aligner("Local",synml)
        final=tf.concat([conv_final_g,conv_final_l],1)

        with tf.variable_scope("Combine_Renormalize",reuse=tf.AUTO_REUSE):
            bl=tf.concat([bls,blhs],1)
            #bl=bls-blhs
            theta_bl=get_variable_by_shape((21*2,1),"theta_bl")
            theta_bl=tf.matmul(bl,theta_bl)
            x=tf.concat([x,theta_bl],1)
            theta=get_variable_by_shape((4,49*2),"theta")
            bias=get_variable_by_shape((1,49*2),"b")
            theta_2=tf.matmul(x,theta)+bias
            print(theta_2)
            theta_2=tf.reshape(theta_2,(-1,49*2,1))
            theta_2=tf.tile(theta_2,[1,1,dim*2])
            print(theta_2)
            final=final*theta_2
            print(final)

        flat=tf.layers.flatten(final)
        #flat_w=tf.layers.flatten(W)
        #flat=tf.concat([flat,flat_w],1)
        print(flat)
        dense=tf.layers.dense(flat,2048,kernel_regularizer=reg_l2,bias_regularizer=reg_l2)
        dense_2=tf.layers.dense(dense,1024,kernel_regularizer=reg_l2,bias_regularizer=reg_l2)
        dense_3=tf.layers.dense(dense_2,512,kernel_regularizer=reg_l2,bias_regularizer=reg_l2)
        
        logits_pred=tf.layers.dense(dense_3,4,name="Predictions")
        entropy=tf.nn.sparse_softmax_cross_entropy_with_logits(logits=logits_pred,labels=y)
        print(entropy)
        
        #weights = tf.trainable_variables() # all vars of your graph
        #regl1 = tf.contrib.layers.apply_regularization(reg_l1, weights)
        reg_losses = tf.get_collection(tf.GraphKeys.REGULARIZATION_LOSSES)
        reg_constant =0.00000001
        loss=tf.reduce_mean(entropy)+reg_constant * sum(reg_losses)
        #loss=tf.reduce_mean(entropy)
        optimizer=tf.train.RMSPropOptimizer(lr)
        #optimizer=tf.train.AdamOptimizer()
        t_op=optimizer.minimize(loss)
        acc=tf.math.in_top_k(tf.cast(logits_pred,tf.float32),y,1)
        accuracy=tf.reduce_mean(tf.cast(acc,tf.float32))
        init=tf.global_variables_initializer()
        saver=tf.train.Saver()

    for node in (synmg,synml,bls,blhs,gl,dps,dphs,dis,lr,y):
        g.add_to_collection("input_nodes",node)

    for node in (loss,t_op,accuracy,init):
        g.add_to_collection("output_nodes",node)

    return g,saver
