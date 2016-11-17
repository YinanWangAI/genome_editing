"""
Deep rank: Rank sgRNAs by Convolutional Neural Network

Now we only use sequence as input, we could consider GC content, position on
the protein, structure etc.
"""

import numpy as np
import scipy.stats
import tensorflow as tf


# Input
def generate_input(seqs, score):
    dataset = []
    for i, seq in enumerate(seqs):
        encoding_seq = one_hot_encoding(seq)
        dataset.append([encoding_seq, score[i]])
    return np.array(dataset)


def one_hot_encoding(seq):
    encoding_dict = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    seq = seq.upper()
    encoding_mat = np.zeros(shape=(4, len(seq)))
    for i, bp in enumerate(seq):
        encoding_mat[encoding_dict[bp]][i] = 1
    return encoding_mat


def split_data_by_gene(dataset, test_gene):
    train_data = []
    test_data = []
    for record in dataset:
        if record[-1] == test_gene:
            test_data.append(record)
        else:
            train_data.append(record)
    train_data = np.array(train_data)
    test_data = np.array(test_data)
    return train_data, test_data


def split_data_random(dataset, train_ratio=0.6, valid_ratio=0.2,
                      test_ratio=0.2):
    assert train_ratio + valid_ratio + test_ratio == 1, 'Wrong ratio'
    total_num = dataset.shape[0]
    train_num = int(total_num * train_ratio)
    valid_num = int(total_num * valid_ratio)

    train_index = np.random.choice(total_num, train_num, replace=False)
    test_index = np.array(list(set(range(total_num)) - set(train_index)))
    train_data = dataset[train_index]
    valid_index = test_index[:valid_num]
    test_index = test_index[valid_num:]
    valid_data = dataset[valid_index]
    test_data = dataset[test_index]

    return train_data, valid_data, test_data


def permute(x, y):
    assert len(x) == len(y), "The number of features and responses is different"
    permute_index = np.random.choice(len(x), len(x), replace=False)
    return x[permute_index], y[permute_index]


# Inference
def inference(input_tensor, keep_prob):
    with tf.name_scope('conv1'):
        kernel = weight_variable([4, 4, 1, 32])
        conv1 = tf.nn.conv2d(input_tensor, kernel, strides=[1, 1, 1, 1],
                             padding='SAME')
        biases = bias_variable([32])
        conv1_act = tf.nn.relu(conv1 + biases)

    with tf.name_scope('conv2'):
        kernel = weight_variable([4, 2, 32, 64])
        conv2 = tf.nn.conv2d(conv1_act, kernel, strides=[1, 1, 1, 1],
                             padding='SAME')
        biases = bias_variable([64])
        conv2_act = tf.nn.relu(conv2 + biases)

    conv2_dim = int(conv2_act.get_shape()[1] * conv2_act.get_shape()[2] *
                    conv2_act.get_shape()[3])
    conv2_flat = tf.reshape(conv2_act, [-1, conv2_dim])

    hidden1 = fc_layer(conv2_flat, conv2_dim, 1024, layer_name='layer1')

    with tf.name_scope('dropout'):
        dropped = tf.nn.dropout(hidden1, keep_prob)

    y = fc_layer(dropped, 1024, 1, layer_name='readout', act=tf.identity)
    y_hat = 1 / (1 + tf.exp(-y))
    return y_hat


def weight_variable(shape):
    return tf.Variable(tf.truncated_normal(shape=shape, stddev=0.1))


def bias_variable(shape):
    return tf.Variable(tf.constant(0.1, shape=shape))


def fc_layer(input_tensor, input_dim, output_dim, layer_name, act=tf.nn.relu):
    with tf.name_scope(layer_name):
        with tf.name_scope('weights'):
            w = weight_variable([input_dim, output_dim])
        with tf.name_scope('biases'):
            b = bias_variable(shape=[output_dim])
        with tf.name_scope('pre_activate'):
            preact = tf.matmul(input_tensor, w) + b
        activations = act(preact)
    return activations


def loss(y_hat, input_y):
    # y_hat = inference(input_x, keep_prob)
    cov = tf.reduce_sum(tf.matmul(tf.transpose(y_hat - tf.reduce_mean(y_hat)),
                                  input_y - tf.reduce_mean(input_y)))
    sd_y = tf.sqrt(
        tf.matmul(tf.transpose(input_y - tf.reduce_mean(input_y)),
                  input_y - tf.reduce_mean(input_y)))
    sd_y_hat = tf.sqrt(tf.matmul(tf.transpose(y_hat - tf.reduce_mean(y_hat)),
                                 y_hat - tf.reduce_mean(y_hat)))
    pearson_r = cov / (sd_y * sd_y_hat)
    cost_function = -pearson_r
    return cost_function


def train_step(cost_function, lr=1e-4,
               optimizer=tf.train.AdamOptimizer):
    step = optimizer(lr).minimize(cost_function)
    return step


def deep_rank(train_x, train_y, valid_x, valid_y, max_epoch=20, batch_size=100,
              model_save_path='deep_rank_model.ckpt'):
    with tf.Graph().as_default():
        # train_x and train_y
        x = tf.placeholder(tf.float32, [None, 4, 30, 1])
        y = tf.placeholder(tf.float32, [None, 1])
        keep_prob = tf.placeholder(tf.float32)

        # inference model
        y_hat = inference(x, keep_prob)

        # loss
        cost_function = loss(y_hat, y)

        # train_op
        train_op = train_step(cost_function)

        # init
        init = tf.initialize_all_variables()

        # saver
        saver = tf.train.Saver(tf.all_variables())

        # sess
        sess = tf.Session()
        sess.run(init)

        # train parameters
        batch_num = int(len(train_x) / batch_size)

        # training
        for epoch in range(max_epoch):
            train_x, train_y = permute(train_x, train_y)
            for i in range(batch_num):
                batch_x = train_x[(i * batch_size): ((i + 1) * batch_size)]
                batch_y = train_y[(i * batch_size): ((i + 1) * batch_size)]
                feed_dict = {x: batch_x, y: batch_y, keep_prob: 0.5}
                sess.run(train_op, feed_dict=feed_dict)
            print(sess.run(cost_function,
                           feed_dict={x: valid_x, y: valid_y, keep_prob: 1})[0])

        save_path = saver.save(sess, model_save_path)
        print('Save model in {}'.format(save_path))
        sess.close()


# Prediction and evaluation
def predict(model_save_path, input_x, input_y):
    with tf.Graph().as_default():
        x = tf.placeholder(tf.float32, [None, 4, 30, 1])
        y = tf.placeholder(tf.float32, [None, 1])
        keep_prob = tf.placeholder(tf.float32)
        y_hat = inference(x, keep_prob)
        cost_function = loss(y_hat, y)
        saver = tf.train.Saver(tf.all_variables())
        with tf.Session() as sess:
            saver.restore(sess, model_save_path)
            y_pred = sess.run(y_hat, feed_dict={x: input_x, keep_prob: 1})
            pearsonr = sess.run(cost_function,
                                feed_dict={x: input_x, y: input_y,
                                           keep_prob: 1})[0]
        return y_pred, pearsonr


def evaluate(y_true, y_pred):
    return scipy.stats.pearsonr(y_true, y_pred)[0], \
           scipy.stats.spearmanr(y_true, y_pred)[0]
