import os
import itertools
import numpy as np
import torch
from torch import nn
from .utilities import ReaderWriter


class RNNHyperparamConfig:
    def __init__(self, 
                 embed_dim, 
                 z_dim,
                 num_hidden_layers, 
                 bidirection, 
                 p_dropout,     
                 rnn_class,
                 nonlin_func, 
                 l2_reg, 
                 batch_size,
                num_epochs):

        self.embed_dim = embed_dim
        self.z_dim = z_dim
        self.num_hidden_layers = num_hidden_layers
        self.bidirection = bidirection
        self.p_dropout = p_dropout
        self.rnn_class = rnn_class
        self.nonlin_func = nonlin_func
        self.l2_reg = l2_reg
        self.batch_size = batch_size
        self.num_epochs = num_epochs


    def __repr__(self):
        desc = f" embed_dim:{self.embed_dim}\n z_dim:{self.z_dim}\n" \
               f" num_hidden_layers:{self.num_hidden_layers}\n " \
              f"  bidirection:{self.bidirection}\n " \
               f" p_dropout:{self.p_dropout} \n rnn_class:{self.rnn_class} \n nonlin_func:{self.nonlin_func} \n " \
               f" l2_reg:{self.l2_reg} \n batch_size:{self.batch_size} \n num_epochs: {self.num_epochs}"
        return desc

def generate_hyperparam_space(model_name):
    if(model_name == 'PE_RNN'):
        embed_dim = [16,32,64,128]
        z_dim = [16,32,64,128]
        num_hidden_layers = [1, 2, 3]
        bidirection = [True, False]
        p_dropout = [0.1, 0.25, 0.45]
        rnn_class = [nn.GRU, nn.LSTM]
        nonlin_func = [nn.ReLU(), nn.ELU()]
        l2_reg = [1e-4, 1e-5, 1e-6]
        batch_size = [1500, 2500]
        num_epochs = [25]
        opt_lst = [embed_dim, z_dim, 
                   num_hidden_layers, bidirection,
                   p_dropout, rnn_class,
                   nonlin_func, l2_reg, batch_size, num_epochs]
    hyperparam_space = list(itertools.product(*opt_lst))

    return hyperparam_space

def compute_numtrials(prob_interval_truemax, prob_estim):
    """ computes number of trials needed for random hyperparameter search
        see `algorithms for hyperparameter optimization paper
        <https://papers.nips.cc/paper/4443-algorithms-for-hyper-parameter-optimization.pdf>`__
        Args:
            prob_interval_truemax: float, probability interval of the true optimal hyperparam,
                i.e. within 5% expressed as .05
            prob_estim: float, probability/confidence level, i.e. 95% expressed as .95
    """
    n = np.log(1-prob_estim)/np.log(1-prob_interval_truemax)
    return(int(np.ceil(n))+1)


def get_hyperparam_options(prob_interval_truemax, prob_estim, model_name, random_seed=42):
    np.random.seed(random_seed)
    num_trials = compute_numtrials(prob_interval_truemax, prob_estim)
    hyperparam_space = generate_hyperparam_space(model_name)
    if(num_trials > len(hyperparam_space)):
        num_trials = len(hyperparam_space)
    indxs = np.random.choice(len(hyperparam_space), size=num_trials, replace=False)

    if(model_name == 'PE_RNN'):
        hyperconfig_class = RNNHyperparamConfig
    return [hyperconfig_class(*hyperparam_space[indx]) for indx in indxs]



def get_saved_config(config_dir):
    options = ReaderWriter.read_data(os.path.join(config_dir, 'exp_options.pkl'))
    mconfig = ReaderWriter.read_data(os.path.join(config_dir, 'mconfig.pkl'))
    return mconfig, options


def get_index_argmax(score_matrix, target_indx):
    argmax_indx = np.argmax(score_matrix, axis=0)[target_indx]
    return argmax_indx

def get_random_run(num_runs, random_seed=42):
    """Get for each experiment the run number to use for identifying optimal hyperparams
    """
    np.random.seed(random_seed)
#     return np.random.randint(num_runs)
    return 0


def get_best_config_from_hyperparamsearch(hyperparam_search_dir, num_trials=60, num_metrics=5, metric_indx=3):
    """Read best models config from all models tested in hyperparamsearch phase
    Args:
        hyperparam_search_dir: string, path root directory where hyperparam models are stored
        num_trials: int, number of tested models (default 60 based on 0.05 interval and 0.95 confidence interval)
                    see :func: `compute_numtrials`
        metric_indx:int, (default 3) using AUPR as performance metric when num_metrics = 5
                         (default 1) using Spearman correlation coefficient as performance metric when num_metrics = 3
    """
    # determine best config from hyperparam search
#     run_num = get_random_run(num_runs, random_seed=random_seed)
    run_dir = os.path.join(hyperparam_search_dir)

    scores = np.ones((num_trials, num_metrics))*-1
    exist_flag = False

    for config_num in range(num_trials):

        score_file = os.path.join(run_dir, 'config_{}'.format(config_num), 'score_validation.pkl')
        if(os.path.isfile(score_file)):
            try:
                mscore = ReaderWriter.read_data(score_file)
                print(mscore)
                if num_metrics == 5:
                    scores[config_num, 0] = mscore.best_epoch_indx
                    scores[config_num, 1] = mscore.binary_f1
                    scores[config_num, 2] = mscore.macro_f1
                    scores[config_num, 3] = mscore.aupr
                    scores[config_num, 4] = mscore.auc
                    exist_flag = True
                elif num_metrics == 3:
                    scores[config_num, 0] = mscore.best_epoch_indx
                    scores[config_num, 1] = mscore.spearman_corr
                    scores[config_num, 2] = mscore.pearson_corr
                    exist_flag = True
            except Exception as e:
                print(f'exception occured at config_{config_num}')
                continue
        else:
            print("WARNING: hyperparam search dir does not exist: {}".format(score_file))
    if(exist_flag):
        argmax_indx = get_index_argmax(scores, metric_indx)
        mconfig, options = get_saved_config(os.path.join(run_dir, 'config_{}'.format(argmax_indx), 'config'))
        return mconfig, options, argmax_indx, scores
    
    return None

