#!/data/www/crispor/venv/bin/python3
import os, sys
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import numpy as np

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_WEIGHT_DIR = os.path.join(_SCRIPT_DIR, "ABE_Proportion_Weight")

length = 26
window_start = 5
window_size = 8


class Deep_xCas9(object):
    def __init__(self, filter_size, filter_num, node_1=80, node_2=60, l_rate=0.005, window_size=5):
        self.inputs = tf.placeholder(tf.float32, [None, 1, length, 4])
        self.targets = tf.placeholder(tf.float32, [None, 2**window_size - 1])
        self.wow = tf.placeholder(tf.float32, [None, 2**window_size - 1])
        self.possible_labels = tf.placeholder(tf.float32, [None, 2**window_size - 1])
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

        with tf.variable_scope("Output_Layer"):
            W_out = tf.get_variable("W_out", shape=[node_1, 2**window_size - 1])
            B_out = tf.get_variable("B_out", shape=[2**window_size - 1])
            self.outputs = tf.nn.bias_add(tf.matmul(L_fcl1_drop, W_out), B_out)

        self.possible_outputs = tf.nn.softmax(self.outputs)
        self.obj_loss = tf.reduce_mean(
            -tf.reduce_sum(
                self.targets * tf.log(self.possible_outputs)
                - self.targets * tf.log(self.targets),
                reduction_indices=[1],
            )
        )
        self.obj_loss1 = tf.reduce_mean(
            -tf.reduce_sum(
                self.targets * tf.log(self.wow) - self.targets * tf.log(self.targets),
                reduction_indices=[1],
            )
        )
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


def inference_index(orig_seq, req_seq):
    index = []
    for seq_index in range(len(list(req_seq))):
        labels_index = -1
        for ind in range(window_size):
            if req_seq[seq_index][window_start + ind] == orig_seq[seq_index][window_start + ind]:
                pass
            else:
                labels_index += 2 ** (window_size - 1 - ind)
        if labels_index < 0:
            labels_index = 0
        index.append(labels_index)
    return index


def req_seq_produce(seq):
    req_seq = []
    full_seq = []
    for indiv_seq in seq:
        tmp_seq = [indiv_seq]
        for index in range(window_size):
            if indiv_seq[window_start + index] == "A":
                tmp = []
                for tmp_indiv_seq in tmp_seq:
                    tmp.append(
                        tmp_indiv_seq[: window_start + index]
                        + "G"
                        + tmp_indiv_seq[window_start + index + 1 :]
                    )
                    full_seq.append(indiv_seq)
                for sequence in tmp:
                    tmp_seq.append(sequence)
        for req_sequence in tmp_seq:
            if req_sequence != tmp_seq[0]:
                req_seq.append(req_sequence)
    return full_seq, req_seq


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
        model = Deep_xCas9(filter_size, filter_num, node_1, node_2, l_rate, window_size)
        saver = tf.train.Saver()
        saver.restore(sess, os.path.join(weight_dir, checkpoint))

    return {"sess": sess, "model": model, "graph": graph}


def predict(handle, sequences, batch_size=1024):
    """Predict edit-outcome proportions for a list of 26-bp sequences.

    For each input sequence, enumerates every variant produced by converting any
    subset of A bases to G within positions [window_start, window_start+window_size).
    Returns a list parallel to `sequences`; each entry is a list of (variant_seq, probability)
    tuples. Probabilities are normalized across the entire predict() call (sums to 1.0
    over the full batch), matching the original script's behavior — call with one
    sequence at a time for per-input normalization.
    """
    full_seq, req_seq = req_seq_produce(list(sequences))
    if len(full_seq) == 0:
        return [[] for _ in sequences]

    X = preprocess_seq(full_seq)
    labels = inference_index(full_seq, req_seq)

    sess = handle["sess"]
    model = handle["model"]
    graph = handle["graph"]

    TEST_Z = np.zeros((X.shape[0], 2**window_size - 1), dtype=float)
    with graph.as_default():
        for i in range(int(np.ceil(X.shape[0] / float(batch_size)))):
            chunk = X[i * batch_size : (i + 1) * batch_size]
            TEST_Z[i * batch_size : (i + 1) * batch_size] = sess.run(
                model.possible_outputs,
                feed_dict={model.inputs: chunk, model.is_training: False},
            )

    total = sum(TEST_Z[i][labels[i]] for i in range(len(TEST_Z)))
    if total > 0:
        TEST_Z = TEST_Z / total

    results = [[] for _ in sequences]
    v_idx = 0
    for input_idx, indiv_seq in enumerate(sequences):
        a_count = sum(1 for j in range(window_size) if indiv_seq[window_start + j] == "A")
        n_variants = (2 ** a_count - 1) if a_count > 0 else 0
        for _ in range(n_variants):
            results[input_idx].append((req_seq[v_idx], float(TEST_Z[v_idx][labels[v_idx]])))
            v_idx += 1
    return results


def close_model(handle):
    handle["sess"].close()
