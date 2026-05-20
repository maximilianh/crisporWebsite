#!/data/www/crispor/venv/bin/python3
import os, sys
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import numpy as np

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_WEIGHT_DIR = os.path.join(_SCRIPT_DIR, "CBE_Efficiency_Weight")

length = 24


class Deep_xCas9(object):
    def __init__(self, filter_size, filter_num, node_1=80, node_2=60, l_rate=0.005):
        self.inputs = tf.placeholder(tf.float32, [None, 1, length, 4])
        self.targets = tf.placeholder(tf.float32, [None, 1])
        self.is_training = tf.placeholder(tf.bool)

        def create_new_conv_layer(
            input_data, num_input_channels, num_filters, filter_shape, pool_shape, name
        ):
            conv_filt_shape = [
                filter_shape[0],
                filter_shape[1],
                num_input_channels,
                num_filters,
            ]

            weights = tf.Variable(
                tf.truncated_normal(conv_filt_shape, stddev=0.03), name=name + "_W"
            )
            bias = tf.Variable(tf.truncated_normal([num_filters]), name=name + "_b")

            out_layer = tf.nn.conv2d(input_data, weights, [1, 1, 1, 1], padding="VALID")
            out_layer += bias
            out_layer = tf.layers.dropout(tf.nn.relu(out_layer), 0.3, self.is_training)
            return out_layer

        L_filter_num = 4
        L_inputs = self.inputs
        L_pool_0 = create_new_conv_layer(
            L_inputs,
            L_filter_num,
            filter_num[0] * 3,
            [1, filter_size[0]],
            [1, 2],
            name="conv1",
        )

        with tf.variable_scope("Fully_Connected_Layer1"):
            layer_node_0 = int((length - filter_size[0]) / 1) + 1
            node_num_0 = layer_node_0 * filter_num[0] * 3

            L_flatten_0 = tf.reshape(L_pool_0, [-1, node_num_0])
            W_fcl1 = tf.get_variable("W_fcl1", shape=[node_num_0, node_1])
            B_fcl1 = tf.get_variable("B_fcl1", shape=[node_1])
            L_fcl1_pre = tf.nn.bias_add(tf.matmul(L_flatten_0, W_fcl1), B_fcl1)
            L_fcl1 = tf.nn.relu(L_fcl1_pre)
            L_fcl1_drop = tf.layers.dropout(L_fcl1, 0.3, self.is_training)

        with tf.variable_scope("Fully_Connected_Layer2"):
            W_fcl2 = tf.get_variable("W_fcl2", shape=[node_1, node_2])
            B_fcl2 = tf.get_variable("B_fcl2", shape=[node_2])
            L_fcl2_pre = tf.nn.bias_add(tf.matmul(L_fcl1_drop, W_fcl2), B_fcl2)
            L_fcl2 = tf.nn.relu(L_fcl2_pre)
            L_fcl2_drop = tf.layers.dropout(L_fcl2, 0.3, self.is_training)

        with tf.variable_scope("Output_Layer"):
            W_out = tf.get_variable("W_out", shape=[node_2, 1])
            B_out = tf.get_variable("B_out", shape=[1])
            self.outputs = tf.nn.bias_add(tf.matmul(L_fcl2_drop, W_out), B_out)

        self.obj_loss = tf.reduce_mean(tf.square(self.targets - self.outputs))
        self.optimizer = tf.train.AdamOptimizer(l_rate).minimize(self.obj_loss)


def preprocess_seq(data):
    DATA_X = np.zeros((len(data), 1, length, 4), dtype=int)
    for l in range(len(data)):
        for i in range(length):
            c = data[l][i]
            if c in "Aa":
                DATA_X[l, 0, i, 0] = 1
            elif c in "Cc":
                DATA_X[l, 0, i, 1] = 1
            elif c in "Gg":
                DATA_X[l, 0, i, 2] = 1
            elif c in "Tt":
                DATA_X[l, 0, i, 3] = 1
            else:
                raise ValueError("Non-ATGC character %r at position %d in %r" % (c, i, data[l]))
    return DATA_X


def _find_checkpoint(weight_dir):
    for name in os.listdir(weight_dir):
        if name.endswith(".meta"):
            return name[:-5]
    raise FileNotFoundError("No .meta checkpoint found in %s" % weight_dir)


def load_model(weight_dir=_DEFAULT_WEIGHT_DIR):
    """Build the graph, restore weights, and return a handle for repeated predict() calls."""
    checkpoint = _find_checkpoint(weight_dir)
    parts = checkpoint.split("-")
    parsed = []
    for v in parts:
        try:
            parsed.append(int(v))
        except ValueError:
            try:
                parsed.append(float(v))
            except ValueError:
                parsed.append(v)

    (
        filter_size_1, filter_size_2, filter_size_3,
        filter_num_1, filter_num_2, filter_num_3,
        l_rate, load_episode, node_1, node_2,
    ) = parsed[2:]

    filter_size = [filter_size_1, filter_size_2, filter_size_3]
    filter_num = [filter_num_1, filter_num_2, filter_num_3]

    conf = tf.ConfigProto()
    conf.gpu_options.allow_growth = True
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", "0")

    graph = tf.Graph()
    with graph.as_default():
        sess = tf.Session(config=conf)
        model = Deep_xCas9(filter_size, filter_num, node_1, node_2, l_rate)
        saver = tf.train.Saver()
        saver.restore(sess, os.path.join(weight_dir, checkpoint))

    return {"sess": sess, "model": model, "graph": graph}


def predict(handle, sequences, batch_size=500):
    """Run inference on a list of 24-bp sequences. Returns a 1D numpy array of efficiency scores (percent)."""
    X = preprocess_seq(sequences)
    sess = handle["sess"]
    model = handle["model"]
    graph = handle["graph"]

    out = np.zeros((X.shape[0],), dtype=float)
    with graph.as_default():
        for i in range(int(np.ceil(X.shape[0] / float(batch_size)))):
            chunk = X[i * batch_size : (i + 1) * batch_size]
            preds = sess.run(model.outputs, feed_dict={model.inputs: chunk, model.is_training: False})
            out[i * batch_size : (i + 1) * batch_size] = 100.0 * preds[:, 0]
    return out


def close_model(handle):
    handle["sess"].close()
