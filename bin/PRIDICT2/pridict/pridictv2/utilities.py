import os
import shutil
import pickle
import string
import torch
import numpy as np
import scipy
from scipy.stats import gaussian_kde
import pandas as pd
from matplotlib import pyplot as plt

class ContModelScore:
    def __init__(self, best_epoch_indx, spearman_corr, pearson_corr):
        self.best_epoch_indx = best_epoch_indx
        self.spearman_corr = spearman_corr
        self.pearson_corr = pearson_corr
        self.spearman_lst = []
        self.pearson_lst = []

    def __repr__(self):
        desc = " best_epoch_indx:{}\n spearman_corr:{}\n pearson_corr:{}\n" \
               "".format(self.best_epoch_indx, self.spearman_corr, self.pearson_corr)
        return desc
def get_performance_multidata_results(target_dir, num_runs, dsettype, ref_run, dsetnames, 
                            outcome_name='averageedited',task_type='continuous', switch=False):

    if task_type == 'categ':
        metric_names = ('auc', 'aupr', 'macro_f1')
    elif task_type == 'continuous':
        metric_names = ('spearman_corr', 'pearson_corr')
    elif task_type == 'ordinal':
        metric_names = ('mae', 'mse')

    num_metrics = len(metric_names)

    outcome_name_indices_map ={'averageedited':0,
                               'averageunedited':1,
                               'averageindel':2}
    outcome_indx = outcome_name_indices_map[outcome_name]
    res_lst = []
    if dsettype in {'train', 'validation'} and ref_run is None:
        prefix = 'train_val'
    else:
        prefix = 'test'
    for i_data, dsetname in enumerate(dsetnames):
        all_perf = {}
        perf_dict = [{} for i in range(num_metrics)]
        for run_num in range(num_runs):
            
            if ref_run is not None:
                runname = 'run_{}_{}'.format(ref_run, run_num)

            else:
                runname = 'run_{}'.format(run_num)
            if switch:
                run_dir = os.path.join(target_dir, 
                                    runname,
                                    '{}'.format(prefix))
            else:
                run_dir = os.path.join(target_dir,
                                    '{}'.format(prefix),
                                    runname)

            score_file = os.path.join(run_dir, 'score_{}.pkl'.format(dsettype))
            if os.path.isfile(score_file):
                mscore = ReaderWriter.read_data(score_file)

                if task_type == 'continuous' and switch == False: # case of neural models
                    dsetmscore = mscore.modelscores_lst[i_data]
                    perf_dict[0][runname] = dsetmscore.spearman_lst[outcome_indx]
                    perf_dict[1][runname] = dsetmscore.pearson_lst[outcome_indx]

        perf_df_lst = []
        for i in range(num_metrics):
            all_perf = perf_dict[i]
            all_perf_df = pd.DataFrame(all_perf, index=[f'{metric_names[i]}_{dsetname}_{outcome_name}'])
            median = all_perf_df.median(axis=1)
            mean = all_perf_df.mean(axis=1)
            stddev = all_perf_df.std(axis=1)
            all_perf_df['mean'] = mean
            all_perf_df['median'] = median
            all_perf_df['stddev'] = stddev
            perf_df_lst.append(all_perf_df.sort_values('mean', ascending=False))
        dset_perf_df = pd.concat(perf_df_lst, axis=0)
        res_lst.append(dset_perf_df)
    return res_lst

def build_performance_multidata_dfs(target_dir, num_runs, dsettype, outcome_name, task_type, dset_names, ref_run=None, switch=False):
    target_dir = create_directory(target_dir)
    return get_performance_multidata_results(target_dir, num_runs, dsettype, ref_run, dset_names, outcome_name=outcome_name, task_type=task_type, switch=switch)

def get_performance_results(target_dir, num_runs, dsettype, ref_run, outcome_name='averageedited',task_type='continuous', switch=False):

    if task_type == 'categ':
        metric_names = ('auc', 'aupr', 'macro_f1')
    elif task_type == 'continuous':
        metric_names = ('spearman_corr', 'pearson_corr')
    elif task_type == 'ordinal':
        metric_names = ('mae', 'mse')

    num_metrics = len(metric_names)
    all_perf = {}
    perf_dict = [{} for i in range(num_metrics)]
    outcome_name_indices_map ={'averageedited':0,
                               'averageunedited':1,
                               'averageindel':2}
    outcome_indx = outcome_name_indices_map[outcome_name]

    if dsettype in {'train', 'validation'} and ref_run is None:
        prefix = 'train_val'
    else:
        prefix = 'test'

    for run_num in range(num_runs):
        
        if ref_run is not None:
            runname = 'run_{}_{}'.format(ref_run, run_num)

        else:
            runname = 'run_{}'.format(run_num)
        if switch:
            run_dir = os.path.join(target_dir, 
                                   runname,
                                   '{}'.format(prefix))
        else:
            run_dir = os.path.join(target_dir,
                                '{}'.format(prefix),
                                runname)

        score_file = os.path.join(run_dir, 'score_{}.pkl'.format(dsettype))
        print(score_file)
        if os.path.isfile(score_file):
            mscore = ReaderWriter.read_data(score_file)


            if task_type == 'categ':
                perf_dict[0][runname] = mscore.auc
                perf_dict[1][runname] = mscore.aupr
                perf_dict[2][runname] = mscore.macro_f1

            elif task_type == 'continuous' and switch == False: # case of neural models
                perf_dict[0][runname] = mscore.spearman_lst[outcome_indx]
                perf_dict[1][runname] = mscore.pearson_lst[outcome_indx]
            elif task_type == 'continuous' and switch == True: # case of baselines
                perf_dict[0][runname] = mscore.spearman_corr
                perf_dict[1][runname] = mscore.pearson_corr

            elif task_type == 'ordinal':
                perf_dict[0][runname] = mscore.mae
                perf_dict[1][runname] = mscore.mse

    perf_df_lst = []
    for i in range(num_metrics):
        all_perf = perf_dict[i]
        all_perf_df = pd.DataFrame(all_perf, index=[f'{metric_names[i]}_{outcome_name}'])
        median = all_perf_df.median(axis=1)
        mean = all_perf_df.mean(axis=1)
        stddev = all_perf_df.std(axis=1)
        all_perf_df['mean'] = mean
        all_perf_df['median'] = median
        all_perf_df['stddev'] = stddev
        perf_df_lst.append(all_perf_df.sort_values('mean', ascending=False))
    
    return pd.concat(perf_df_lst, axis=0)


def build_performance_dfs(target_dir, num_runs, dsettype, outcome_name, task_type, ref_run=None, switch=False):
    target_dir = create_directory(target_dir)
    return get_performance_results(target_dir, num_runs, dsettype, ref_run, outcome_name=outcome_name, task_type=task_type, switch=switch)

def update_Adamoptimizer_lr_momentum_(optm, lr, momen):
    """in-place update for learning rate and momentum for Adam optimizer"""
    for pg in optm.param_groups:
        pg['lr'] = lr
        pg['betas'] = (momen, pg['betas'][-1])

def compute_lr_scheduler(l0, lmax, num_iter, annealing_percent):
    num_annealing_iter = np.floor(annealing_percent * num_iter)
    num_iter_upd = num_iter - num_annealing_iter
    
    x = [0, np.ceil(num_iter_upd/2.0), num_iter_upd, num_iter]
    y = [l0, lmax, l0, l0/100.0]
    tck = scipy.interpolate.splrep(x, y, k=1, s=0)
    
    xnew = np.arange(0, num_iter)
    lrates = scipy.interpolate.splev(xnew, tck, der=0)
    return lrates

def compute_momentum_scheduler(momen_0, momen_max, num_iter, annealing_percent):
    num_annealing_iter = np.floor(annealing_percent * num_iter)
    num_iter_upd = num_iter - num_annealing_iter
    
    x = [0, np.ceil(num_iter_upd/2.0), num_iter_upd, num_iter]
    y = [momen_max, momen_0, momen_max, momen_max]
    tck = scipy.interpolate.splrep(x, y, k=1, s=0)
    
    xnew = np.arange(0, num_iter)
    momentum_vals = scipy.interpolate.splev(xnew, tck, der=0)
    return momentum_vals


def build_predictions_df(seq_ids, true_score, pred_score, y_ref_names, dset_names=None):
         
    
    seqid_inpseq_df = pd.DataFrame(seq_ids)
    seqid_inpseq_df.columns = ['seq_id']

    if dset_names is not None:
        seqid_inpseq_df['dataset_name'] = dset_names
    
    
    target_names = ['averageedited', 'averageunedited', 'averageindel']

    assert (len(y_ref_names) > 0 and len(y_ref_names) <= 3), f'# of target outcomes should be > 0 and not exceed 3!. Possible outcome names are:\n {target_names}'

    df_dict = {}


    if true_score is not None:
        true_score_arr = np.array(true_score)
        if len(true_score_arr.shape) == 1: # (nsamples,)
            true_score_arr = true_score_arr.reshape(-1,1)
        
        num_targets = true_score_arr.shape[-1]
        assert num_targets <= 3, 'number of targets should not exceed three outcomes: averageedited, averageunedtied, averageidnel!'
        true_scores_dict = {}
        for i in range (num_targets):
            target_name = y_ref_names[i]
            true_scores_dict[f'true_{target_name}'] = true_score_arr[:, i]
        df_dict.update(true_scores_dict)

    pred_score_arr = np.array(pred_score)
    if len(pred_score_arr.shape) == 1: # (nsamples,)
        pred_score_arr = pred_score_arr.reshape(-1, 1)
    
    num_targets = pred_score_arr.shape[-1]
    assert num_targets in {1,2, 3}, '# of predicted outcomes should be at max 3'
    pred_scores_dict = {}
    if num_targets > 1:
        for i in range (num_targets):
            target_name = target_names[i]
            pred_scores_dict[f'pred_{target_name}'] = pred_score_arr[:, i]
    elif num_targets == 1:
        target_name = y_ref_names[0]
        pred_scores_dict[f'pred_{target_name}'] = pred_score_arr[:, 0]
 

    # print('true_score_arr.shape:', true_score_arr.shape)
    # print('pred_score_arr.shape:',pred_score_arr.shape)
    df_dict.update(pred_scores_dict)
    predictions_df = pd.concat([seqid_inpseq_df, pd.DataFrame(df_dict)], axis=1)

    # print('pred.shape:', predictions_df.shape)
    if dset_names is not None:
        predictions_df.drop_duplicates(['seq_id', 'dataset_name'], inplace=True)
        # print('after removing duplicates -> pred.shape:', predictions_df.shape)

    return predictions_df

def dump_dict_content(dsettype_content_map, dsettypes, desc, wrk_dir):
    for dsettype in dsettypes:
        path = os.path.join(wrk_dir, '{}_{}.pkl'.format(desc, dsettype))
        ReaderWriter.dump_data(dsettype_content_map[dsettype], path)

class ReaderWriter(object):
    """class for dumping, reading and logging data"""
    def __init__(self):
        pass

    @staticmethod
    def dump_data(data, file_name, mode="wb"):
        """dump data by pickling
           Args:
               data: data to be pickled
               file_name: file path where data will be dumped
               mode: specify writing options i.e. binary or unicode
        """
        with open(file_name, mode) as f:
            pickle.dump(data, f)

    @staticmethod
    def read_data(file_name, mode="rb"):
        """read dumped/pickled data
           Args:
               file_name: file path where data will be dumped
               mode: specify writing options i.e. binary or unicode
        """
        with open(file_name, mode) as f:
            data = pickle.load(f)
        return(data)

    @staticmethod
    def dump_tensor(data, file_name):
        """
        Dump a tensor using PyTorch's custom serialization. Enables re-loading the tensor on a specific gpu later.
        Args:
            data: Tensor
            file_name: file path where data will be dumped
        Returns:
        """
        torch.save(data, file_name)

    @staticmethod
    def read_tensor(file_name, device):
        """read dumped/pickled data
           Args:
               file_name: file path where data will be dumped
               device: the gpu to load the tensor on to
        """
        data = torch.load(file_name, map_location=device)
        return data

    @staticmethod
    def write_log(line, outfile, mode="a"):
        """write data to a file
           Args:
               line: string representing data to be written out
               outfile: file path where data will be written/logged
               mode: specify writing options i.e. append, write
        """
        with open(outfile, mode) as f:
            f.write(line)

    @staticmethod
    def read_log(file_name, mode="r"):
        """write data to a file
           Args:
               line: string representing data to be written out
               outfile: file path where data will be written/logged
               mode: specify writing options i.e. append, write
        """
        with open(file_name, mode) as f:
            for line in f:
                yield line

def switch_layer_to_traineval_mode(model, target_layer, activate_train=True):
    """
    Target a layer and switch to train or eval mode
    """
    for child_name, child in model.named_children():
        if isinstance(child, target_layer):
            print(child_name, '=>', target_layer, 'is switching training mode to ', activate_train)
            if activate_train:
                child.train()
            else:
                child.eval()
            
        else:
            switch_layer_to_traineval_mode(child, target_layer, activate_train=activate_train)

def grad_track_hook(tensor_name):
    def print_hook(grad):
        pass
        # print('grad for ', tensor_name, ' is computed with grad_shape:', grad.shape, ' and grad_nrom:', grad.norm())
        # print('*'*15)
    return print_hook

def require_nonleaf_grad(v, tensor_name):
    v.retain_grad()
    v.register_hook(grad_track_hook(tensor_name))

##### utilities to retune trained model on another smaller set of experiments ####

def get_trainable_params(models):
    models_param = []
    for layer, layer_name in models:
        for params in layer.parameters():
            if params.requires_grad:
                models_param.extend(list(params))
    return models_param

def get_num_trainable_params(models):
    trainable_params_count = 0.
    for layer, layer_name in models:
        for param_name, params in layer.named_parameters():
            if params.requires_grad:
                num_params = params.numel()
                print(layer_name, param_name, num_params)
                trainable_params_count+= num_params
    return trainable_params_count

def freeze_layers(models, trainable_layernames):
    for layer, layer_name in models:
        for param_name, params in layer.named_parameters():
            if layer_name not in trainable_layernames:
                params.requires_grad = False
            print(layer_name, param_name, params.requires_grad)
        print()

def create_directory(folder_name, directory="current"):
    """create directory/folder (if it does not exist) and returns the path of the directory
       Args:
           folder_name: string representing the name of the folder to be created
       Keyword Arguments:
           directory: string representing the directory where to create the folder
                      if `current` then the folder will be created in the current directory
    """
    if directory == "current":
        path_current_dir = os.path.dirname(__file__)  # __file__ refers to utilities.py
    else:
        path_current_dir = directory
    path_new_dir = os.path.join(path_current_dir, folder_name)
    if not os.path.exists(path_new_dir):
        os.makedirs(path_new_dir)
    return(path_new_dir)


def get_device(to_gpu, index=0):
    is_cuda = torch.cuda.is_available()
    if(is_cuda and to_gpu):
        target_device = 'cuda:{}'.format(index)
    else:
        target_device = 'cpu'
    return torch.device(target_device)


def report_available_cuda_devices():
    if(torch.cuda.is_available()):
        n_gpu = torch.cuda.device_count()
        print('number of GPUs available:', n_gpu)
        for i in range(n_gpu):
            print("cuda:{}, name:{}".format(i, torch.cuda.get_device_name(i)))
            device = torch.device('cuda', i)
            get_cuda_device_stats(device)
            print()
    else:
        print("no GPU devices available!!")

def get_cuda_device_stats(device):
    print('total memory available:', torch.cuda.get_device_properties(device).total_memory/(1024**3), 'GB')
    print('total memory allocated on device:', torch.cuda.memory_allocated(device)/(1024**3), 'GB')
    print('max memory allocated on device:', torch.cuda.max_memory_allocated(device)/(1024**3), 'GB')
    print('total memory cached on device:', torch.cuda.memory_cached(device)/(1024**3), 'GB')
    print('max memory cached  on device:', torch.cuda.max_memory_cached(device)/(1024**3), 'GB')

def compute_harmonic_mean(a, b):
    if a==0 and b==0:
        return 0.
    return 2*a*b/(a+b)
    
def compute_spearman_corr(pred_score, ref_score):
    # return scipy.stats.kendalltau(pred_score, ref_score)
    return scipy.stats.spearmanr(pred_score, ref_score)

def compute_pearson_corr(pred_score, ref_score):
    # return scipy.stats.kendalltau(pred_score, ref_score)
    return scipy.stats.pearsonr(pred_score, ref_score)

def restrict_grad_(mparams, mode, limit):
    """clamp/clip a gradient in-place
    """
    if(mode == 'clip_norm'):
        __, maxl = limit
        torch.nn.utils.clip_grad_norm_(mparams, maxl, norm_type=2) # l2 norm clipping
    elif(mode == 'clamp'): # case of clamping
        minl, maxl = limit
        for param in mparams:
            if param.grad is not None:
                param.grad.data.clamp_(minl, maxl)
def check_na(df):
    assert df.isna().any().sum() == 0

def perfmetric_report_cont(pred_score, ref_score, epoch_loss, epoch, outlog):
    """
    Args:
        pred_score: list, (nsamples,) or (nsamples, num_targets)
        ref_score: list, (nsamples,) or (nsamples, num_targets)
    """
    lsep = "\n"
    report = "Epoch: {}".format(epoch) + lsep
  
    pred_score_arr = np.array(pred_score)
    ref_score_arr = np.array(ref_score)
    # print('pred_score_arr.shape:', pred_score_arr.shape)
    # print('ref_score_arr.shape:', ref_score_arr.shape)
    
    # (nsamples, distribution)
    if len(pred_score_arr.shape) > 1: # case of predicting distribution
        target_names = ['averageedited', 'averageunedited', 'averageindel']
        num_targets = pred_score_arr.shape[-1] # number of classes
        spear_lst = []
        pear_lst = []
        for target_indx in range(num_targets):
            spearman_corr, pvalue_spc = compute_spearman_corr(pred_score_arr[:,target_indx], ref_score_arr[:,target_indx])
            pearson_corr, pvalue_prc = compute_pearson_corr(pred_score_arr[:,target_indx], ref_score_arr[:,target_indx])
            spear_lst.append(spearman_corr)
            pear_lst.append(pearson_corr)
            target_name = target_names[target_indx]
            report += f"Spearman correlation score {target_name}:{spearman_corr}   pvalue:{pvalue_spc}" + lsep
            report += f"Pearson correlation score {target_name}:{pearson_corr}   pvalue:{pvalue_prc}" + lsep    
            report += "-"*15 + lsep
        # compute overall correlation
        spearman_corr = np.mean(spear_lst[0:1]+spear_lst[-1:])
        pearson_corr = np.mean(pear_lst[0:1]+pear_lst[-1:])
        modelscore = ContModelScore(epoch, spearman_corr, pearson_corr)
        # TODO: clean this hack
        modelscore.spearman_lst = spear_lst
        modelscore.pearson_lst = pear_lst
    # (nsamples,)
    else:
        spearman_corr, pvalue_spc = compute_spearman_corr(pred_score_arr, ref_score_arr)
        pearson_corr, pvalue_prc = compute_pearson_corr(pred_score_arr, ref_score_arr)
        report += f"Spearman correlation score:{spearman_corr}    pvalue:{pvalue_spc}" + lsep
        report += f"Pearson correlation score:{pearson_corr}    pvalue:{pvalue_prc}" + lsep    
        modelscore = ContModelScore(epoch, spearman_corr, pearson_corr)
        # TODO: clean this hack
        modelscore.spearman_lst = [spearman_corr]
        modelscore.pearson_lst = [pearson_corr]

    report += f"epoch average batch loss:{epoch_loss}" + lsep
    report += "-"*15 + lsep
    ReaderWriter.write_log(report, outlog)
    return modelscore

def compute_performance_multidata_from_df(df, dataset_names, target_names=['averageedited', 'averageunedited', 'averageindel']):
    # (nsamples, distribution)
    df = df.copy()
    # make sure there are no duplicates 
    print('df.shape:', df.shape)
    df.drop_duplicates(['seq_id', 'dataset_name'], inplace=True)
    print('after deduplicating, df.shape:', df.shape)
    modelscores_lst = []

    lsep='\n'
    report=''
    for dname in dataset_names:
        cond = df['dataset_name'] == dname
        spear_lst = []
        pear_lst = []
        report += f"peformance of dataset_name:{dname}" + lsep
        for target_indx in range(len(target_names)):
            target_name = target_names[target_indx]
            pred_score_arr = df.loc[cond, f'pred_{target_name}'].values
            ref_score_arr = df.loc[cond, f'true_{target_name}'].values
            spearman_corr, pvalue_spc = compute_spearman_corr(pred_score_arr, ref_score_arr)
            pearson_corr, pvalue_prc = compute_pearson_corr(pred_score_arr, ref_score_arr)
            spear_lst.append(spearman_corr)
            pear_lst.append(pearson_corr)
            report += f"Spearman correlation score {target_name}:{spearman_corr}   pvalue:{pvalue_spc}" + lsep
            report += f"Pearson correlation score {target_name}:{pearson_corr}   pvalue:{pvalue_prc}" + lsep    
            report += "-"*15 + lsep
        # compute overall correlation
        # spearman_corr = np.mean(spear_lst[0:1]+spear_lst[-1:])
        # pearson_corr = np.mean(pear_lst[0:1]+pear_lst[-1:])
        spearman_corr = np.mean(spear_lst[0:1])
        pearson_corr = np.mean(pear_lst[0:1])
        modelscore = ContModelScore(-1, spearman_corr, pearson_corr)
        # TODO: clean this hack
        modelscore.spearman_lst = spear_lst
        modelscore.pearson_lst = pear_lst
        modelscores_lst.append(modelscore)
    # create avarege modelscore across the datasets
    mean_spearmancorr = np.mean([modelscore.spearman_corr for modelscore in modelscores_lst])
    mean_pearsoncorr = np.mean([modelscore.pearson_corr for modelscore in modelscores_lst])
    avg_modelscore = ContModelScore(-1, mean_spearmancorr, mean_pearsoncorr)
    avg_modelscore.dataset_names = dataset_names
    avg_modelscore.modelscores_lst = modelscores_lst

    return avg_modelscore, report

def perfmetric_report_multidata_cont(pred_score, ref_score, dset_ids, seq_ids, epoch_loss, epoch, outlog):
    """
    Args:
        pred_score: list, (nsamples,) or (nsamples, num_targets)
        ref_score: list, (nsamples,) or (nsamples, num_targets)
    """
    lsep = "\n"
    report = "Epoch: {}".format(epoch) + lsep
  
    pred_score_arr = np.array(pred_score)
    ref_score_arr = np.array(ref_score)

    # print('pred_score_arr.shape:', pred_score_arr.shape)
    # print('ref_score_arr.shape:', ref_score_arr.shape)
    df = pd.DataFrame({'dataset_name':dset_ids,
                       'seq_id':seq_ids})
    print('__perfmetric_report_multidata_cont__')
    print('# of rows:', df.shape[0])

    dataset_names = df['dataset_name'].unique().tolist()

    # (nsamples, distribution)
    # TODO: just send the target_names as variable to the function as we are doing in :func:`compute_performance_multidata_from_df`
    if pred_score_arr.shape[-1]  == 3: # case of predicting distribution
        target_names = ['averageedited', 'averageunedited', 'averageindel']
    elif pred_score_arr.shape[-1]  == 2:
        target_names = ['averageedited', 'averageunedited']
    else:
        target_names = ['averageedited']

    # print('target_names:', target_names)

    for target_index, target_name in enumerate(target_names):
        df[f'pred_{target_name}'] = pred_score_arr[:,target_index]
        df[f'ref_{target_name}'] = ref_score_arr[:,target_index]

    df.drop_duplicates(['seq_id', 'dataset_name'], inplace=True)
    print('updated # of rows:', df.shape[0])

    modelscores_lst = []

    for dname in dataset_names:
        cond = df['dataset_name'] == dname
        spear_lst = []
        pear_lst = []
        report += f"peformance of dataset_name:{dname}" + lsep
        for target_indx in range(len(target_names)):
            target_name = target_names[target_indx]
            pred_score_arr = df.loc[cond, f'pred_{target_name}'].values
            ref_score_arr = df.loc[cond, f'ref_{target_name}'].values
            spearman_corr, pvalue_spc = compute_spearman_corr(pred_score_arr, ref_score_arr)
            pearson_corr, pvalue_prc = compute_pearson_corr(pred_score_arr, ref_score_arr)
            spear_lst.append(spearman_corr)
            pear_lst.append(pearson_corr)
            report += f"Spearman correlation score {target_name}:{spearman_corr}   pvalue:{pvalue_spc}" + lsep
            report += f"Pearson correlation score {target_name}:{pearson_corr}   pvalue:{pvalue_prc}" + lsep    
            report += "-"*15 + lsep
        # compute overall correlation
        # spearman_corr = np.mean(spear_lst[0:1]+spear_lst[-1:])
        # pearson_corr = np.mean(pear_lst[0:1]+pear_lst[-1:])
        spearman_corr = np.mean(spear_lst[0:1])
        pearson_corr = np.mean(pear_lst[0:1])
        modelscore = ContModelScore(epoch, spearman_corr, pearson_corr)
        # TODO: clean this hack
        modelscore.spearman_lst = spear_lst
        modelscore.pearson_lst = pear_lst
        modelscores_lst.append(modelscore)
    # create avarege modelscore across the datasets
    mean_spearmancorr = np.mean([modelscore.spearman_corr for modelscore in modelscores_lst])
    mean_pearsoncorr = np.mean([modelscore.pearson_corr for modelscore in modelscores_lst])
    avg_modelscore = ContModelScore(epoch, mean_spearmancorr, mean_pearsoncorr)
    avg_modelscore.dataset_names = dataset_names
    avg_modelscore.modelscores_lst = modelscores_lst

    report += f"epoch average batch loss:{epoch_loss}" + lsep
    report += "-"*15 + lsep
    ReaderWriter.write_log(report, outlog)
    return avg_modelscore

def plot_loss(epoch_loss_avgbatch, wrk_dir):
    dsettypes = epoch_loss_avgbatch.keys()
    for dsettype in dsettypes:
        plt.figure(figsize=(9, 6))
        plt.plot(epoch_loss_avgbatch[dsettype], 'r')
        plt.xlabel("number of epochs")
        plt.ylabel("negative loglikelihood cost")
        plt.legend(['epoch batch average loss'])
        plt.savefig(os.path.join(wrk_dir, os.path.join(dsettype + ".pdf")))
        plt.close()


def plot_xy(x, y, xlabel, ylabel, legend, fname, wrk_dir):
    plt.figure(figsize=(9, 6))
    plt.plot(x, y, 'r')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if legend:
        plt.legend([legend])
    plt.savefig(os.path.join(wrk_dir, os.path.join(fname + ".pdf")))
    plt.close()

def delete_directory(directory):
    if(os.path.isdir(directory)):
        shutil.rmtree(directory)
        
def transform_genseq_upper(df, columns):
    for colname in columns:
        df[colname] = df[colname].str.upper()
    return df

def plot_corr(y_pred, y_ref, model_name, dsettype, sort_pts_density=False, fig_size=(5,5),figdir=None):
    if sort_pts_density:
        xy = np.vstack([y_pred,y_ref])
        z = gaussian_kde(xy)(xy)
        # Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        corr_x, corr_y, corr_z = y_pred[idx], y_ref[idx], z[idx]
    else:
        corr_z = 'blue'
    R = compute_spearman_corr(y_pred,y_ref)[0]
    r = compute_pearson_corr(y_pred, y_ref)[0]
    fig, ax = plt.subplots(1,1,figsize = fig_size,dpi=100)
    plt.rcParams['axes.linewidth'] = 1.5
    ax.scatter(y_ref,y_pred, c=corr_z, s=2, label = f'{model_name}:{dsettype} set', alpha=0.3)
#     ax.set_title(f'{model_name} - Correlation\n{dsettype} Data', fontsize=16)
    ax.set_ylabel('Predicted efficiency (%)', fontsize=14)
    ax.set_xlabel('Measured efficiency (%)', fontsize=14)
    m, b = np.polyfit(y_ref,y_pred, 1)
    ax.plot(y_ref, m*y_ref + b, color='black', alpha=0.6)
    ax.set_ylim(0,)
    ax.text(0.7,0.15, 'R = '+str(R)[:4]+'\nr = '+str(r)[:4], transform=ax.transAxes, fontsize=14,
            verticalalignment='top')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='best', fontsize=10)
    if figdir:
        figpth = os.path.join(figdir, f"corr_{model_name}_{dsettype}.png")
        fig.savefig(os.path.join(figpth),bbox_inches='tight')
        plt.close()

