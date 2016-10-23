import numpy as np
from sklearn.preprocessing import OneHotEncoder
import tensorflow as tf


def encode_dna(seq, n_sample, n_feature,
               encode_dic={'A': 0, 'T': 1, 'C': 2, 'G': 3}):
    seq_int = []
    for i in range(n_sample):
        seq_int.append(char2int(seq[i], encode_dic))
    seq_int = np.asarray(seq_int)

    enc = OneHotEncoder(sparse=False, n_values=4)
    data_onehot = enc.fit_transform(seq_int)

    out = [np.reshape(x, (n_feature, 4)) for x in data_onehot]
    return np.asarray(out)


def char2int(s, encode_dic):
    out = []
    for char in s:
        out.append(encode_dic[char])
    return np.asarray(out)