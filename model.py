dim=10
n=3
fl=2*n+1
maxbl=29

import tensorflow as tf

def create_model():
    tf.compat.v1.reset_default_graph()
    g=tf.Graph()
    with g.as_default():
        synmg=tf.compat.v1.placeholder(dtype=tf.float64,shape=(None,2*n+1,2*n+1,2),name="Synteny_matrix_placeholder_Global")
        synml=tf.compat.v1.placeholder(dtype=tf.float64,shape=(None,2*n+1,2*n+1,2),name="Synteny_matrix_placeholder_Local")
        pfam=tf.compat.v1.placeholder(dtype=tf.float64,shape=(None,2*n+1,2*n+1,1),name="Pfam_matrix_placeholder")
        bls=tf.compat.v1.placeholder(dtype=tf.float64,shape=(None,maxbl),name="Species_Branch_Length_Placeholder")
        blhs=tf.compat.v1.placeholder(dtype=tf.float64,shape=(None,maxbl),name="Homology_Species_Branch_Length_Placeholder")
        #gl=tf.placeholder(dtype=tf.float64,shape=(None,1),name="Mean_gene_length")
        dps=tf.compat.v1.placeholder(dtype=tf.float64,shape=(None,1),name="mca_species_distance")
        dphs=tf.compat.v1.placeholder(dtype=tf.float64,shape=(None,1),name="mca_homology_species_distance")
        dis=tf.compat.v1.placeholder(dtype=tf.float64,shape=(None,1),name="total_distance")
        lr=tf.compat.v1.placeholder(dtype=tf.float64,shape=(),name="learning_rate")
        y=tf.compat.v1.placeholder(dtype=tf.int32,shape=(None),name="labels")
        
        lrs=tf.compat.v1.summary.scalar("Learning_Rate",lr)
        
        x=tf.concat([dps,dps-dphs,dis],1,name="Create_train_vector")

        print(synmg,"\n",synml,"\n",bls,"\n",blhs,"\n",dps,"\n",dphs,"\n",dis,"\n",x)

        reg_l2=tf.keras.regularizers.l2(0.5 * (0.001))
        reg_l1 = tf.keras.regularizers.l1(l=0.005)

        def get_variable_by_shape(shape,name):
            f=tf.compat.v1.get_variable(name,shape=shape,initializer=tf.compat.v1.glorot_uniform_initializer(),dtype=tf.float64,regularizer=reg_l2)
            return f

        def create_synteny_aligner(name,synm):
            with tf.compat.v1.variable_scope(name+"Synteny_Aligner",reuse=tf.compat.v1.AUTO_REUSE):

                fconv=get_variable_by_shape((2,2,2,dim),"fconv")
                conv=tf.nn.conv2d(input=synm,filters=fconv,strides=(1,1,1,1),padding="VALID",name="Conv_aligner")

                fconv_1=get_variable_by_shape((2,2,dim,dim*2),"fconv_1")
                conv_1=tf.nn.conv2d(input=conv,filters=fconv_1,strides=(1,1,1,1),padding="VALID",name="Conv_aligner_1")
                
                fxconv=get_variable_by_shape((fl,2,dim*2),"fxconv")
                x_conv=tf.reshape(synm,(-1,fl*fl,2))
                x_conv=tf.nn.conv1d(input=x_conv,filters=fxconv,stride=fl,padding="SAME",name="row_aligner")
                
                y_conv=tf.reshape(tf.transpose(a=synm,perm=(0,2,1,3)),(-1,fl*fl,2))
                fyconv=get_variable_by_shape((fl,2,dim*2),"fyconv")
                y_conv=tf.nn.conv1d(input=y_conv,filters=fyconv,stride=fl,padding="SAME",name="column_aligner")
                
                
                wconv=get_variable_by_shape((fl,fl,2,dim*2*10),"wconv")
                w_conv=tf.nn.conv2d(input=synm,filters=wconv,strides=(1,1,1,1),padding="VALID",name="Global_Aligner_1")
                
                conv_1=tf.reshape(conv_1,(-1,25,dim*2))
                x_conv=tf.reshape(x_conv,(-1,fl,dim*2))
                y_conv=tf.reshape(y_conv,(-1,fl,dim*2))
                w_conv=tf.reshape(w_conv,(-1,10,dim*2))
                conv_final=tf.concat([conv_1,x_conv,y_conv,w_conv],1,name="Concatenate_All_Alignments")
            return conv_final

        with tf.compat.v1.variable_scope("Pfam",reuse=tf.compat.v1.AUTO_REUSE):
            wconv_pfam=get_variable_by_shape((fl,fl,1,dim*2*10),"wconv_pfam")
            w_conv_pfam=tf.nn.conv2d(input=pfam,filters=wconv_pfam,strides=(1,1,1,1),padding="VALID",name="Global_Aligner_pfam")
            w_conv_pfam=tf.reshape(w_conv_pfam,(-1,10,dim*2))

        conv_final_g=create_synteny_aligner("Global_",synmg)
        conv_final_l=create_synteny_aligner("Local_",synml)
        final=tf.concat([conv_final_g,conv_final_l,w_conv_pfam],1)
        #final=conv_final_l

        with tf.compat.v1.variable_scope("Combine_Renormalize",reuse=tf.compat.v1.AUTO_REUSE):
            bl=tf.concat([bls,blhs],1)
            #bl=bls-blhs
            theta_bl=get_variable_by_shape((maxbl*2,1),"theta_bl")
            theta_bl=tf.matmul(bl,theta_bl)
            x=tf.concat([x,theta_bl],1)
            theta=get_variable_by_shape((4,108),"theta")
            bias=get_variable_by_shape((1,108),"b")
            theta_2=tf.matmul(x,theta)+bias
            theta_2=tf.reshape(theta_2,(-1,108,1))
            theta_2=tf.tile(theta_2,[1,1,dim*2])
            final=final*theta_2
            print(final)

        flat=tf.compat.v1.layers.flatten(final)

        zero=tf.constant(0.0,dtype=tf.float64)
        diff=dps-dphs
        print(diff)
        diff_2=tf.cast(tf.equal(diff,zero),tf.float64)
        print(diff_2)
        diff=tf.tile(diff_2,[1,dim*10])

        flat=tf.concat([flat,diff],1)
        print(flat)
        #dense=tf.layers.dense(flat,2048,kernel_regularizer=reg_l2,bias_regularizer=reg_l2)
        #dense_2=tf.layers.dense(dense,1024,kernel_regularizer=reg_l2,bias_regularizer=reg_l2)
        dense_3=tf.compat.v1.layers.dense(flat,512,kernel_regularizer=reg_l2,bias_regularizer=reg_l2)
        
        logits_pred=tf.compat.v1.layers.dense(dense_3,3,name="Predictions")
        print(logits_pred)
        entropy=tf.nn.sparse_softmax_cross_entropy_with_logits(logits=logits_pred,labels=y)
        print(entropy)
        
        #weights = tf.trainable_variables() # all vars of your graph
        #regl1 = tf.contrib.layers.apply_regularization(reg_l1, weights)
        reg_losses = tf.compat.v1.get_collection(tf.compat.v1.GraphKeys.REGULARIZATION_LOSSES)
        reg_constant =0.00000001
        loss=tf.reduce_mean(input_tensor=entropy)+reg_constant * sum(reg_losses)
        #loss=tf.reduce_mean(entropy)
        optimizer=tf.compat.v1.train.RMSPropOptimizer(lr)
        #optimizer=tf.train.AdamOptimizer()
        losses=tf.compat.v1.summary.scalar("Loss",loss)
        t_op=optimizer.minimize(loss)
        acc=tf.math.in_top_k(predictions=tf.cast(logits_pred,tf.float32),targets=y,k=1)
        accuracy=tf.reduce_mean(input_tensor=tf.cast(acc,tf.float32))
        accs=tf.compat.v1.summary.scalar("Accuracy",accuracy)
        summary=tf.compat.v1.summary.merge_all()
        init=tf.compat.v1.global_variables_initializer()
        saver=tf.compat.v1.train.Saver()

    for node in (synmg,synml,pfam,bls,blhs,dps,dphs,dis,lr,y):
        g.add_to_collection("input_nodes",node)

    for node in (loss,t_op,accuracy,init,summary):
        g.add_to_collection("output_nodes",node)

    return g,saver
