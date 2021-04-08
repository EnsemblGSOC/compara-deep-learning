import tensorflow as tf
dim = 10
n = 3
fl = 2*n+1
maxbl = 29


def create_model():
    tf.reset_default_graph()
    g = tf.Graph()
    with g.as_default():
        synmg = tf.placeholder(dtype=tf.float64, shape=(
            None, 2*n+1, 2*n+1, 2), name="Synteny_matrix_placeholder_Global")
        synml = tf.placeholder(dtype=tf.float64, shape=(
            None, 2*n+1, 2*n+1, 2), name="Synteny_matrix_placeholder_Local")
        pfam = tf.placeholder(dtype=tf.float64, shape=(
            None, 2*n+1, 2*n+1, 1), name="Pfam_matrix_placeholder")
        bls = tf.placeholder(dtype=tf.float64, shape=(
            None, maxbl), name="Species_Branch_Length_Placeholder")
        blhs = tf.placeholder(dtype=tf.float64, shape=(
            None, maxbl), name="Homology_Species_Branch_Length_Placeholder")
        # gl=tf.placeholder(dtype=tf.float64,shape=(None,1),name="Mean_gene_length")
        dps = tf.placeholder(dtype=tf.float64, shape=(
            None, 1), name="mca_species_distance")
        dphs = tf.placeholder(dtype=tf.float64, shape=(
            None, 1), name="mca_homology_species_distance")
        dis = tf.placeholder(dtype=tf.float64, shape=(
            None, 1), name="total_distance")
        lr = tf.placeholder(dtype=tf.float64, shape=(), name="learning_rate")
        y = tf.placeholder(dtype=tf.int32, shape=(None), name="labels")

        # lrs = tf.summary.scalar("Learning_Rate", lr)

        x = tf.concat([dps, dps-dphs, dis], 1, name="Create_train_vector")

        print(synmg, "\n", synml, "\n", bls, "\n", blhs,
              "\n", dps, "\n", dphs, "\n", dis, "\n", x)

        reg_l2 = tf.contrib.layers.l2_regularizer(0.001)
        # reg_l1 = tf.contrib.layers.l1_regularizer(scale=0.005, scope=None)

        def get_variable_by_shape(shape, name):
            f = tf.get_variable(name, shape=shape,
                                initializer=tf.glorot_uniform_initializer(),
                                dtype=tf.float64,
                                regularizer=reg_l2)
            return f

        def create_synteny_aligner(name, synm):
            with tf.variable_scope(
                                  name+"Synteny_Aligner", reuse=tf.AUTO_REUSE):

                fconv = get_variable_by_shape((2, 2, 2, dim), "fconv")
                conv = tf.nn.conv2d(synm, fconv, (1, 1, 1, 1),
                                    padding="VALID", name="Conv_aligner")

                fconv_1 = get_variable_by_shape((2, 2, dim, dim*2), "fconv_1")
                conv_1 = tf.nn.conv2d(
                    conv, fconv_1, (1, 1, 1, 1), padding="VALID",
                    name="Conv_aligner_1")

                fxconv = get_variable_by_shape((fl, 2, dim*2), "fxconv")
                x_conv = tf.reshape(synm, (-1, fl*fl, 2))
                x_conv = tf.nn.conv1d(
                    x_conv, fxconv, stride=fl, padding="SAME",
                    name="row_aligner")

                y_conv = tf.reshape(tf.transpose(
                    synm, (0, 2, 1, 3)), (-1, fl*fl, 2))
                fyconv = get_variable_by_shape((fl, 2, dim*2), "fyconv")
                y_conv = tf.nn.conv1d(
                    y_conv, fyconv, stride=fl, padding="SAME",
                    name="column_aligner")

                wconv = get_variable_by_shape((fl, fl, 2, dim*2*10), "wconv")
                w_conv = tf.nn.conv2d(
                    synm, wconv, (1, 1, 1, 1), padding="VALID",
                    name="Global_Aligner_1")

                conv_1 = tf.reshape(conv_1, (-1, 25, dim*2))
                x_conv = tf.reshape(x_conv, (-1, fl, dim*2))
                y_conv = tf.reshape(y_conv, (-1, fl, dim*2))
                w_conv = tf.reshape(w_conv, (-1, 10, dim*2))
                conv_final = tf.concat(
                    [conv_1, x_conv, y_conv, w_conv], 1,
                    name="Concatenate_All_Alignments")
            return conv_final

        with tf.variable_scope("Pfam", reuse=tf.AUTO_REUSE):
            wconv_pfam = get_variable_by_shape(
                (fl, fl, 1, dim*2*10), "wconv_pfam")
            w_conv_pfam = tf.nn.conv2d(
                pfam, wconv_pfam, (1, 1, 1, 1), padding="VALID",
                name="Global_Aligner_pfam")
            w_conv_pfam = tf.reshape(w_conv_pfam, (-1, 10, dim*2))

        conv_final_g = create_synteny_aligner("Global_", synmg)
        conv_final_l = create_synteny_aligner("Local_", synml)
        final = tf.concat([conv_final_g, conv_final_l, w_conv_pfam], 1)
        # final=conv_final_l

        with tf.variable_scope("Combine_Renormalize", reuse=tf.AUTO_REUSE):
            bl = tf.concat([bls, blhs], 1)
            # bl=bls-blhs
            theta_bl = get_variable_by_shape((maxbl*2, 1), "theta_bl")
            theta_bl = tf.matmul(bl, theta_bl)
            x = tf.concat([x, theta_bl], 1)
            theta = get_variable_by_shape((4, 108), "theta")
            bias = get_variable_by_shape((1, 108), "b")
            theta_2 = tf.matmul(x, theta)+bias
            theta_2 = tf.reshape(theta_2, (-1, 108, 1))
            theta_2 = tf.tile(theta_2, [1, 1, dim*2])
            final = final*theta_2
            print(final)

        flat = tf.layers.flatten(final)

        zero = tf.constant(0.0, dtype=tf.float64)
        diff = dps-dphs
        print(diff)
        diff_2 = tf.cast(tf.equal(diff, zero), tf.float64)
        print(diff_2)
        diff = tf.tile(diff_2, [1, dim*10])

        flat = tf.concat([flat, diff], 1)
        print(flat)
        # dense=tf.layers.dense(flat,2048,kernel_regularizer=reg_l2,bias_regularizer=reg_l2)
        # dense_2=tf.layers.dense(dense,1024,kernel_regularizer=reg_l2,bias_regularizer=reg_l2)
        dense_3 = tf.layers.dense(
            flat, 512, kernel_regularizer=reg_l2, bias_regularizer=reg_l2)

        logits_pred = tf.layers.dense(dense_3, 3, name="Predictions")
        print(logits_pred)
        entropy = tf.nn.sparse_softmax_cross_entropy_with_logits(
            logits=logits_pred, labels=y)
        print(entropy)

        # weights = tf.trainable_variables() # all vars of your graph
        # regl1 = tf.contrib.layers.apply_regularization(reg_l1, weights)
        reg_losses = tf.get_collection(tf.GraphKeys.REGULARIZATION_LOSSES)
        reg_constant = 0.00000001
        loss = tf.reduce_mean(entropy)+reg_constant * sum(reg_losses)
        # loss=tf.reduce_mean(entropy)
        optimizer = tf.train.RMSPropOptimizer(lr)
        # optimizer=tf.train.AdamOptimizer()
        # losses = tf.summary.scalar("Loss", loss)
        t_op = optimizer.minimize(loss)
        acc = tf.math.in_top_k(tf.cast(logits_pred, tf.float32), y, 1)
        accuracy = tf.reduce_mean(tf.cast(acc, tf.float32))
        # accs = tf.summary.scalar("Accuracy", accuracy)
        summary = tf.summary.merge_all()
        init = tf.global_variables_initializer()
        saver = tf.train.Saver()

    for node in (synmg, synml, pfam, bls, blhs, dps, dphs, dis, lr, y):
        g.add_to_collection("input_nodes", node)

    for node in (loss, t_op, accuracy, init, summary):
        g.add_to_collection("output_nodes", node)

    return g, saver
