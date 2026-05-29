import glob, os
from joblib import load
from sklearn.ensemble import GradientBoostingRegressor
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import numpy as np 
import pandas as pd

PACKAGE_DIRECTORY = os.path.dirname(os.path.realpath(__file__))
MODEL_DIR = PACKAGE_DIRECTORY + '/saved_models/'
try:
    assert os.path.exists(MODEL_DIR)
except:
    raise FileNotFoundError(f'Could not find saved models. Looked in: {MODEL_DIR}')
POS_MODELS = {}
TOTAL_EFFICIENCY_MODELS = {}
POS_FRACS = {
    2: 0.05184948248657137,
    3: 0.19673314654134325,
    4: 0.443757534968822,
    5: 0.9119183639518773,
    6: 1.0,
    7: 0.7846071775331374,
    8: 0.41387793695191166,
    9: 0.14388146489658163,
    10: 0.07095196797688848
}
EDITOR_TARGET_BASE = {
    'BE4': 'C',
    'FNLS': 'C',
    'CBE': 'C',
    'ABE8e': 'A',
    'ABE20m': 'A',
    'ABE': 'A'
}

def load_models():
    global POS_MODELS
    global TOTAL_EFFICIENCY_MODELS
    editors = [os.path.basename(d) for d in glob.glob(MODEL_DIR + '/*')]
    for editor in editors:
        print(f'Loading models for {editor}')
        model_paths = sorted(glob.glob(MODEL_DIR + f'{editor}/pos*'), key=lambda x: int(x.split('.')[0].split('_')[-1]))
        pos_list = []
        pos_models = []
        for path in model_paths:
            pos = int(path.split('.')[0].split('_')[-1])
            model = load(path)
            pos_list.append(pos)
            pos_models.append(model)
        POS_MODELS[editor] = {}
        POS_MODELS[editor]['pos_list'] = pos_list
        POS_MODELS[editor]['pos_models'] = pos_models

        TOTAL_EFFICIENCY_MODELS[editor] = load(glob.glob(MODEL_DIR + f'{editor}/global*')[0])

    print(f'Finished loading models')

# Microhomology features
def generate_microhomology_matrix(seq):
    base_eq_mat = np.zeros((len(seq), len(seq)))
    for i in list(range(len(seq)))[::-1]:
        for j in list(range(len(seq)))[::-1]:
            eq = seq[i] == seq[j]
            if eq:
                base_eq_mat[i, j] = 1
                if (i < len(seq)-1) and (j < len(seq)-1):
                    base_eq_mat[i, j] += base_eq_mat[i+1, j+1]
    cols = [s for s in seq]
    rows = [s for s in seq]
    bbb = pd.DataFrame(base_eq_mat, columns=cols, index=rows)
    return bbb

def find_microhomology_about_cutsite(seq, pam=20, window=20):
    cutsite = pam-3
    mh_mat = generate_microhomology_matrix(seq).iloc[max(0, cutsite-window):cutsite, cutsite:cutsite+window]
    x, y = np.unravel_index(np.argmax(mh_mat.values, axis=None), mh_mat.shape)
    max_mh_len = mh_mat.iloc[x, y]
    max_dist_between_mh = mh_mat.shape[0]-x + y#cutsite + y - x
    return max_mh_len, max_dist_between_mh

all_nucs = ['G', 'A', 'C', 'T']
all_dinucs = [f'{a}{b}' for a in all_nucs for b in all_nucs]
def nucleotide_features(seq):
    feature_labels = []
    features = []
    
    # single nucleotides 
    for pos in range(20):
        nuc = seq[pos]
        for b in all_nucs:
            feature_labels.append(f'{b} at pos {pos+1}')
            features.append(1 if nuc == b else 0)
    
    return features, feature_labels

def featurize_20nt_target(target_seq):
    feature_labels = []
    features = []
        
    # target base counts
    for nuc in all_nucs:
        feature_labels.append(f'{nuc} count')
        features.append(sum([n == nuc for n in target_seq]))

    # nucleotide features
    nuc_features, nuc_labels = nucleotide_features(target_seq)
    features += nuc_features
    feature_labels += nuc_labels

    # melting temperatures
    features.append(mt.Tm_NN(target_seq))
    feature_labels.append('Melting Temperature')

    return feature_labels, features

def scale_zscores(zscores, pos_list, mean, std):
    scaled = []
    for z, pos in zip(zscores, pos_list):
        mean_frac = POS_FRACS[pos]
        pos_mean = mean * mean_frac
        if (z is None) or np.isnan(z):
            s = z 
        else:
            s = (z * std) + pos_mean
            s = np.clip(s, 0, 1)
        scaled.append(s)
    return scaled

def predict(target_seq, mean=None, std=None, editor=None):
    if len(target_seq) != 20:
        raise ValueError(f'Input sequence is {len(target_seq)}, but require 20')
    feature_labels, features = featurize_20nt_target(target_seq)
    features = np.array(features).reshape(1, -1)

    # Get models for editor
    if editor is None:
        editor = 'CBE'
    elif editor not in POS_MODELS:
        raise ValueError(f'Unknown editor: {editor}')
    pos_models = POS_MODELS[editor]['pos_models']
    pos_list = POS_MODELS[editor]['pos_list']

    # Predict per position
    predictions = [model.predict(features)[0] for model in pos_models]
    
    # Rescale if required
    if (mean is not None) and (std is not None):
        predictions = scale_zscores(predictions, pos_list, mean, std)
        
    ret = [(pos, pred) for pos, pred in zip(pos_list, predictions)]
    
    # Overwrite predictions where there is no target base at position
    target_base = EDITOR_TARGET_BASE[editor]
    ret = [(pos, pred) if target_seq[pos-1] == target_base else (pos, None) for pos, pred in zip(pos_list, predictions)]
    return ret

def predict_total(target_seq, mean=None, std=None, editor=None):
    if len(target_seq) != 20:
        raise ValueError(f'Input sequence is {len(target_seq)}, but require 20')
    feature_labels, features = featurize_20nt_target(target_seq)
    features = np.array(features).reshape(1, -1)

    # Get model for editor
    if editor is None:
        editor = 'CBE'
    elif editor not in POS_MODELS:
        raise ValueError(f'Unknown editor: {editor}')
    model = TOTAL_EFFICIENCY_MODELS[editor]

    prediction = model.predict(features)[0]

    # Rescale if required
    if (mean is not None) and (std is not None):
        prediction = scale_zscores([prediction], [6], mean, std)[0]

    return prediction



def predict_batch_fasta(fasta_path, output_path=None, mean=None, std=None, editor=None):
    # Get models for editor
    if editor is None:
        editor = 'CBE'
    elif editor not in POS_MODELS:
        raise ValueError(f'Unknown editor: {editor}')
    pos_models = POS_MODELS[editor]['pos_models']
    pos_list = POS_MODELS[editor]['pos_list']

    # Predict per position for each guide
    zscore_predictions = []
    guides = []
    for record in SeqIO.parse(fasta_path, 'fasta'):
        if (len(record.seq) != 20) or ('N' in record.seq):
            continue
        preds = [pred for pos, pred in predict(record.seq, editor=editor)]
        zscore_predictions.append(preds)
        guides.append(record.id)

    df = pd.DataFrame(zscore_predictions, index=guides, columns=[f'Pos {pos} Zscore' for pos in pos_list])

    # Scale if required
    if (mean is not None) and (std is not None):
        scaled_predictions = [scale_zscores(p, pos_list, mean, std) for p in zscore_predictions]
        scaled_df = pd.DataFrame(scaled_predictions, index=guides, columns=[f'Pos {pos} Scaled' for pos in pos_list])
        df = pd.concat([df, scaled_df], axis=1)
        df['Mean'] = mean
        df['Std Dev'] = std

    df.index.name = 'Oligo Id'

    if output_path is not None:
        df.to_csv(output_path)

    return df


load_models()
