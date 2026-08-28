import os
from itertools import cycle
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit, GroupKFold
import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader
from .utilities import ContModelScore

class PEDataTensor(Dataset):

    def __init__(self, 
                 X_init_nucl, 
                 X_init_proto, 
                 X_init_pbs,
                 X_init_rt,
                 X_mut_nucl,
                 X_mut_pbs,
                 X_mut_rt,
                 seqlevel_feat,
                 seqlevel_feat_colnames,
                 y_score, 
                 x_init_len,
                 x_mut_len,
                 indx_seqid_map):
        # B: batch elements; T: sequence length
        # tensor.float32, (B, T), (sequence characters are mapped to 0-3) and 4 for padded characters
        self.X_init_nucl = X_init_nucl
        self.X_init_proto = X_init_proto
        self.X_init_pbs = X_init_pbs
        self.X_init_rt = X_init_rt

        self.X_mut_nucl = X_mut_nucl
        self.X_mut_pbs = X_mut_pbs
        self.X_mut_rt = X_mut_rt
        
        self.seqlevel_feat = seqlevel_feat
        self.seqlevel_feat_colnames = seqlevel_feat_colnames
        # tensor.float32, (B, outcome_dim)
        self.y_score = y_score  
        # tensor.int32, (B,), (length of each sequence)
        self.x_init_len = x_init_len 
        self.x_mut_len = x_mut_len
        # dictionary {indx:seq_id}
        self.indx_seqid_map = indx_seqid_map
        self.num_samples = self.X_init_nucl.size(0)  # int, number of sequences

    def __getitem__(self, indx):
        if self.y_score is None:
            y_val = -1.
        else:
            y_val = self.y_score[indx]

        # return(self.X_init_nucl[indx],
        #        self.X_init_proto[indx],
        #        self.X_init_pbs[indx],
        #        self.X_init_rt[indx],
        #        self.X_mut_nucl[indx],
        #        self.X_mut_pbs[indx],
        #        self.X_mut_rt[indx],
        #        self.x_init_len[indx],
        #        self.x_mut_len[indx],
        #        y_val,
        #        indx, 
        #        self.indx_seqid_map[indx])
        return(self.X_init_nucl[indx],
               self.X_init_proto[indx],
               self.X_init_pbs[indx],
               self.X_init_rt[indx],
               self.X_mut_nucl[indx],
               self.X_mut_pbs[indx],
               self.X_mut_rt[indx],
               self.x_init_len[indx],
               self.x_mut_len[indx],
               self.seqlevel_feat[indx],
               y_val,
               indx, 
               self.indx_seqid_map[indx])


    def __len__(self):
        return(self.num_samples)

class PartitionDataTensor(Dataset):

    def __init__(self, pe_datatensor, partition_ids, dsettype, run_num):
        self.pe_datatensor = pe_datatensor  # instance of :class:`PEDatatensor`
        self.partition_ids = partition_ids  # list of sequence indices
        self.dsettype = dsettype  # string, dataset type (i.e. train, validation, test)
        self.run_num = run_num  # int, run number
        self.num_samples = len(self.partition_ids[:])  # int, number of sequences in the partition

    def __getitem__(self, indx):
        target_id = self.partition_ids[indx]
        return self.pe_datatensor[target_id]

    def __len__(self):
        return(self.num_samples)

def extend_matrix(mat, ext_mat_shape, fill_val):
    assert len(mat.shape) == 2
    ext_mat = torch.full(ext_mat_shape, fill_val)
    return torch.cat([mat, ext_mat], axis=1)


class MinMaxNormalizer:
    def __init__(self, include_MFE=True, include_addendumfeat=True, normalizer_dict=None):
        self.length_norm = ['Correction_Length', 'Correction_Length_effective',
                            'RToverhangmatches', 'RToverhanglength', 'RTlength','PBSlength']
        self.mfe_norm = ['MFE_protospacer', 
                        'MFE_protospacer_scaffold', 
                        'MFE_extension', 
                        'MFE_extension_scaffold', 
                        'MFE_protospacer_extension_scaffold',
                        'MFE_rt',
                        'MFE_pbs']
        self.mt_norm = ['RTmt', 'RToverhangmt','PBSmt',
                        'protospacermt','extensionmt',
                        'original_base_mt','edited_base_mt']
        
        # to test their contribution
        self.mt_addendum_norm = ['Tm1', 'Tm2', 'Tm2new', 'Tm3', 'Tm4', 'TmD']
        self.basecount_norm = ['nGCcnt1', 'nGCcnt2', 'nGCcnt3', 'fGCcont1', 'fGCcont2', 'fGCcont3']
        
        if normalizer_dict:
            self.normalizer_dict = normalizer_dict

        if include_addendumfeat:
            self.normalizer_info_minmax = [(self.length_norm, 0., 50.), 
                                        (self.mfe_norm, -120., 0.), 
                                        (self.mt_norm, 0., 200.),
                                        (self.mt_addendum_norm, 0., 200.),
                                        (self.basecount_norm, 0., 100.)]
            self.normalizer_info_max = [(self.length_norm, 50.), 
                                        (self.mfe_norm, 120.), 
                                        (self.mt_norm, 200.),
                                        (self.mt_addendum_norm, 200.),
                                        (self.basecount_norm, 100.)]
        else:
            self.normalizer_info_minmax = [(self.length_norm, 0., 50.), 
                                        (self.mfe_norm, -120., 0.), 
                                        (self.mt_norm, 0., 200.)]
            self.normalizer_info_max = [(self.length_norm, 50.), 
                                        (self.mfe_norm, 120.), 
                                        (self.mt_norm, 200.)] 
        if not include_MFE:

            self.normalizer_info_minmax = self.normalizer_info_minmax[:1] + self.normalizer_info_minmax[2:]
            self.normalizer_info_max = self.normalizer_info_max[:1] + self.normalizer_info_max[2:]
            
        self.colnames = self.get_colnames()
  
    def get_colnames(self):
        colnames = []
        for tup in self.normalizer_info_max:
            featnames, __ = tup
            colnames += featnames
        return colnames

    def normalize_cont_cols(self, df, normalize_opt = 'max', suffix=''):
        if normalize_opt == 'max':
            print('--- max normalization ---')
            return self.normalize_cont_cols_max(df, suffix=suffix)
        elif normalize_opt == 'minmax':
            print('--- minmax normalization ---')
            return self.normalize_cont_cols_minmax(df, suffix=suffix)
        elif normalize_opt == 'standardize':
            print('--- standardized normalization ---')
            return self.normalize_cont_cols_meanstd(df, suffix=suffix)
        
    def normalize_cont_cols_meanstd(self, df, suffix=''):
        """inplace min-max normalization of columns"""
        normalizer_info = self.normalizer_info_minmax
        cont_colnames = []
        for colgrp in normalizer_info:
            colnames, __, __ = colgrp
            for colname in colnames:
                m, std = self.normalizer_dict[colname]
                df[colname+suffix] = (df[colname] - m)/std
                cont_colnames.append(colname+suffix)
        return cont_colnames

    def normalize_cont_cols_minmax(self, df, suffix=''):
        """inplace min-max normalization of columns"""
        normalizer_info = self.normalizer_info_minmax
        cont_colnames = []
        for colgrp in normalizer_info:
            colnames, min_val, max_val = colgrp
            for colname in colnames:
                df[colname+suffix] = ((df[colname] - min_val)/(max_val - min_val)).clip(lower=0., upper=1.)
                cont_colnames.append(colname+suffix)
        return cont_colnames

    def normalize_cont_cols_max(self, df, suffix=''):
        """inplace max normalization of columns"""
        normalizer_info = self.normalizer_info_max
        cont_colnames = []
        for colgrp in normalizer_info:
            colnames, max_val = colgrp
            for colname in colnames:
                df[colname+suffix] = df[colname]/max_val
                cont_colnames.append(colname+suffix)
        return cont_colnames

def get_seqlevel_featnames(suffix='_norm'):
    minmax_normalizer = MinMaxNormalizer()
    cont_colnames = minmax_normalizer.get_colnames()
    # by default normalized column names end with _norm suffix
    if suffix is None:
        norm_colnames =  [f'{fname}' for fname in cont_colnames]
    else:
        norm_colnames =  [f'{fname}{suffix}' for fname in cont_colnames]

    seqfeat_cols = [norm_colnames[0]] + \
                   ['Correction_Deletion', 'Correction_Insertion', 'Correction_Replacement'] + \
                   norm_colnames[1:] + ['original_base_mt_nan', 'edited_base_mt_nan']
    return seqfeat_cols

def create_datatensor(data_df, proc_seq_init_df, num_init_cols,  proc_seq_mut_df, num_mut_cols, cont_cols, window=10, y_ref=[]):
    """create a instance of DataTensor from processeed/cleaned dataframe
    
    Args:
        data_df: pandas.DataFrame, dataset
    """
    
    print('--- create_datatensor ---')

    lower_thr = 0
    upper_thr = 150      #    A  C  T  G  -  N
    blank_nucl_token = 5 # B {0, 1, 2, 3, 4, 5}
    blank_pbs_annot_token = 2 # PBS {0,1,2}
    blank_annot_token = 3 # RT_init, RT_mut, Protos {0, 1, 2, 3}

    # process initial sequences
    start_init_seqs = (proc_seq_init_df['start_seq'] - window).clip(lower=lower_thr, upper=None)
    st_init_colindx = start_init_seqs.min()
    
    # end of RT stretch
    end_init_seqs = proc_seq_init_df['end_seq']
    # print('end_init_seqs:\n', end_init_seqs.value_counts())
    # to put assert statment
    assert ((end_init_seqs - start_init_seqs) <= upper_thr).all(), f'Difference between end and start seuqnce should be at most {upper_thr} bp'
    end_init_colindx = end_init_seqs.max()
    upper_init_thr = np.min([st_init_colindx+upper_thr, num_init_cols-1, end_init_colindx + window]) # -1 to compensate for not including end value
    end_init_seqs = (end_init_seqs + window).clip(lower=None, upper=upper_init_thr)
    end_init_colindx = end_init_seqs.max()
    # print('updated end_init_seqs:\n', end_init_seqs.value_counts())
    # print('st_init_colindx:', st_init_colindx)
    # print('end_init_colindx:', end_init_colindx)
    # print('x_init_len_comp:\n', (end_init_seqs - start_init_seqs).value_counts())

    X_init_nucl = torch.from_numpy(proc_seq_init_df[[f'B{i}' for  i in range(st_init_colindx, end_init_colindx+1)]].values).long()
    X_init_proto = torch.from_numpy(proc_seq_init_df[[f'Protos{i}' for  i in range(st_init_colindx, end_init_colindx+1)]].values).long()
    X_init_pbs = torch.from_numpy(proc_seq_init_df[[f'PBS{i}' for  i in range(st_init_colindx, end_init_colindx+1)]].values).long()
    X_init_rt = torch.from_numpy(proc_seq_init_df[[f'RT_init{i}' for  i in range(st_init_colindx, end_init_colindx+1)]].values).long()
    

    # x_init_len = (X_init_nucl != 4).sum(axis=1).long()

    # use ().values for x_init_len and x_mut_len to keep them in array format rather pd.Series
    x_init_len  = (end_init_seqs - start_init_seqs).values
    # print(x_init_len)

    # use blank indicator when nucleotide is N (i.e. blank token)
    ntoken_cond = (X_init_nucl == blank_nucl_token)

    X_init_proto[ntoken_cond] = blank_annot_token
    X_init_rt[ntoken_cond] = blank_annot_token
    X_init_pbs[ntoken_cond] = blank_pbs_annot_token

    # print('X_init_nucl.unique():', X_init_nucl.unique())
    # print('ntoken_cond:', ntoken_cond.unique())
    # print('X_init_proto:', X_init_proto.unique())
    # print('X_init_pbs:', X_init_pbs.unique())
    # print('X_init_rt:', X_init_rt.unique())

    # process mutation sequences
    if 'start_seq' in proc_seq_mut_df:
        start_mut_seqs = (proc_seq_mut_df['start_seq'] - window).clip(lower=lower_thr, upper=None)
    else:
        start_mut_seqs = start_init_seqs
    
    st_mut_colindx = start_mut_seqs.min()

    end_mut_seqs = proc_seq_mut_df['end_seq']
    # print('end_mut_seqs:\n', end_mut_seqs.value_counts())

    assert ((end_mut_seqs - start_mut_seqs) <= upper_thr).all(), f'Difference between end and start seuqnce should be at most {upper_thr} bp'

    end_mut_colindx = end_mut_seqs.max()
    upper_mut_thr = np.min([st_mut_colindx+upper_thr, num_mut_cols-1, end_mut_colindx + window])
    end_mut_seqs = (end_mut_seqs + window).clip(lower=None, upper=upper_mut_thr)
    end_mut_colindx = end_mut_seqs.max()
    # print('updated end_mut_seqs:\n', end_mut_seqs.value_counts())

    # print('st_mut_colindx:', st_mut_colindx)
    # print('end_mut_colindx:', end_mut_colindx)
    # print('x_mut_len_comp:\n', (end_mut_seqs - start_mut_seqs).value_counts())

    X_mut_nucl = torch.from_numpy(proc_seq_mut_df[[f'B{i}' for  i in range(st_mut_colindx, end_mut_colindx+1)]].values).long()
    X_mut_pbs = torch.from_numpy(proc_seq_mut_df[[f'PBS{i}' for  i in range(st_mut_colindx, end_mut_colindx+1)]].values).long()
    X_mut_rt = torch.from_numpy(proc_seq_mut_df[[f'RT_mut{i}' for  i in range(st_mut_colindx, end_mut_colindx+1)]].values).long()

    # x_mut_len = (X_mut_nucl != 4).sum(axis=1).long()
    x_mut_len  = (end_mut_seqs - start_mut_seqs).values

    # use blank indicator when nucleotide is N (i.e. blank token)
    ntoken_cond = (X_mut_nucl == blank_nucl_token)
    X_mut_pbs[ntoken_cond] = blank_pbs_annot_token
    X_mut_rt[ntoken_cond] = blank_annot_token
    # print('X_mut_nucl.unique():', X_mut_nucl.unique())
    # print('ntoken_cond:', ntoken_cond.unique())
    # print('X_mut_pbs:', X_mut_pbs.unique())
    # print('X_mut_rt:', X_mut_rt.unique())

    # harmonize the size of matrices above
    max_num_cols = np.max([end_init_colindx+1, end_mut_colindx+1])
    # print('max_num_cols:', max_num_cols)
        
    #TODO: check if blank_annot_tok or 0 is best
    annot_addendum_tok = blank_annot_token
    # annot_addendum_tok = 0
    if max_num_cols > (end_init_colindx+1):
        num_cols_toadd = max_num_cols - (end_init_colindx + 1)
        bsize = X_init_nucl.shape[0]
        ext_mat_shape = (bsize, num_cols_toadd)
        # mut_cols > init_cols
        X_init_nucl = extend_matrix(X_init_nucl, ext_mat_shape, blank_nucl_token)
        X_init_proto = extend_matrix(X_init_proto, ext_mat_shape, blank_annot_token)
        X_init_pbs = extend_matrix(X_init_pbs, ext_mat_shape, blank_pbs_annot_token)
        X_init_rt = extend_matrix(X_init_rt, ext_mat_shape, blank_annot_token)

    elif max_num_cols > (end_mut_colindx+1):
        num_cols_toadd = max_num_cols - (end_mut_colindx +1)
        bsize = X_mut_nucl.shape[0]
        ext_mat_shape = (bsize, num_cols_toadd)
        # init_cols > mut_cols
        X_mut_nucl = extend_matrix(X_mut_nucl, ext_mat_shape, blank_nucl_token)
        X_mut_pbs = extend_matrix(X_mut_pbs, ext_mat_shape, blank_pbs_annot_token)
        X_mut_rt = extend_matrix(X_mut_rt, ext_mat_shape, blank_annot_token)
    else:
        print('--- aligned initial and mutated sequences ---')

    seq_ids = data_df['seq_id'].values
    indx_seqid_map = {i:seq_ids[i] for i in range(len(seq_ids))}

    if len(y_ref):
        # y_score = torch.from_numpy(data_df['y'].values).reshape(-1,1)
        ycols = [tcol for tcol in y_ref if tcol in data_df]
        # print('ycols:', ycols)
        y_score = torch.from_numpy(data_df[ycols].values)
    else:
        y_score = None

    # get computed features at sequence level
    # seqfeat_cols = ['Correction_Length',
    #                 'Correction_Deletion', 
    #                 'Correction_Insertion', 
    #                 'Correction_Replacement',
    #                 'RToverhangmatches',
    #                 'RToverhanglength',
    #                 'RTlength',
    #                 'PBSlength',
    #                 'MFE_protospacer',
    #                 'MFE_protospacer_scaffold',
    #                 'MFE_extension',
    #                 'MFE_extension_scaffold',
    #                 'MFE_protospacer_extension_scaffold',
    #                 'MFE_rt',
    #                 'MFE_pbs',
    #                 'RTmt',
    #                 'RToverhangmt',
    #                 'PBSmt',
    #                 'protospacermt',
    #                 'extensionmt',
    #                 'original_base_mt',
    #                 'edited_base_mt']
    # seqfeat_cols = ['Correction_Deletion', 
    #                 'Correction_Insertion', 
    #                 'Correction_Replacement']
    # seqfeat_cols += cont_cols
    seqfeat_cols = [cont_cols[0]] + ['Correction_Deletion', 'Correction_Insertion', 'Correction_Replacement'] + cont_cols[1:]

    # case of having indicator variables when melting temperature cannot be computed
    for colname in ('original_base_mt_nan', 'edited_base_mt_nan'):
        if colname in data_df:
            seqfeat_cols.append(colname)

    # Identify the boolean columns
    bool_cols = data_df[seqfeat_cols].select_dtypes(include=['bool']).columns

    # Convert boolean columns to int
    data_df[bool_cols] = data_df[bool_cols].astype(int)
    seqlevel_feat = torch.from_numpy(data_df[seqfeat_cols].values)

    dtensor = PEDataTensor(X_init_nucl, 
                            X_init_proto, 
                            X_init_pbs,
                            X_init_rt,
                            X_mut_nucl,
                            X_mut_pbs,
                            X_mut_rt,
                            seqlevel_feat,
                            seqfeat_cols,
                            y_score, 
                            x_init_len,
                            x_mut_len,
                            indx_seqid_map)
        
    return dtensor


def get_stratified_partitions(group_label, num_folds=5, valid_set_portion=0.1, random_state=42):
    """Generate 5-fold stratified sample of based on passed stratification label y

    Args:
        group_label: group label
    """
    np.random.seed(random_state)
    gkf_trte = GroupKFold(n_splits=num_folds) # split train and test
    gkf_trv = GroupShuffleSplit(n_splits=2, 
                                test_size=valid_set_portion, 
                                random_state=random_state)  # split train and validation
    data_partitions = {}
    X = np.zeros(len(group_label))
    fold_num = 0
    for train_index, test_index in gkf_trte.split(X=X,y=None,groups=group_label):
        
        x_tr = np.zeros(len(train_index))
        y_tr = group_label[train_index]

        for tr_index, val_index in gkf_trv.split(X=x_tr, y=None, groups=y_tr):
            tr_ids = train_index[tr_index]
            val_ids = train_index[val_index]
            data_partitions[fold_num] = {'train': tr_ids,
                                         'validation': val_ids,
                                         'test': test_index}
            
        print("fold_num:", fold_num)
        print('train data')
        print('n=', len(tr_ids))
        # report_label_distrib(y[tr_ids])
        print('validation data')
        print('n=', len(val_ids))
        # report_label_distrib(y[val_ids])
        print('test data')
        print('n=', len(test_index))
        # report_label_distrib(y[test_index])
        print()
        fold_num += 1
        print("-"*25)
    return(data_partitions)

def validate_partitions(data_partitions, sample_ids, valid_set_portion=0.1, test_set_portion=0.2):
    if(not isinstance(sample_ids, set)):
        sample_ids = set(sample_ids)
    num_pts = len(sample_ids)
    print('total number of ids:', num_pts)
    tr_set_accum = set([])
    val_set_accum = set([])
    test_set_accum = set([])
    all_accum = set([])
    for fold_num in data_partitions:
        print('fold_num', fold_num)
        tr_ids = data_partitions[fold_num]['train']
        val_ids = data_partitions[fold_num]['validation']
        te_ids = data_partitions[fold_num]['test']

        tr_val = set(tr_ids).intersection(val_ids)
        tr_te = set(tr_ids).intersection(te_ids)
        te_val = set(te_ids).intersection(val_ids)
        
        tr_size = len(tr_ids) + len(val_ids)
        # assert there is no overlap among train and test partition within a fold
        # print('expected validation set size:', valid_set_portion*tr_size, '; actual validation set size:', len(val_ids))
        # print('expected test set size:', test_set_portion*num_pts, '; actual test set size:', len(te_ids))
        # print()
        # assert there is no interesection between dataset types within a fold
        for s in (tr_val, tr_te, te_val):
            assert len(s) == 0
        s_union = set(tr_ids).union(val_ids).union(te_ids)
        assert len(s_union) == num_pts
        tr_set_accum = tr_set_accum.union(tr_ids)
        val_set_accum = val_set_accum.union(val_ids)
        test_set_accum = test_set_accum.union(te_ids)
        all_accum = all_accum.union(tr_set_accum).union(val_set_accum).union(test_set_accum)
    # verify that assembling test sets from each of the five folds would be equivalent to all drugpair ids
    print('total accum train ids:', len(tr_set_accum))
    print('total accum val ids:', len(val_set_accum))
    print('total accum test ids:', len(test_set_accum))
    print('union of all accum: ', len(all_accum))
    # assert len(test_set_accum) == num_pts
    # assert test_set_accum == sample_ids
    print("passed intersection and overlap test (i.e. train, validation and test sets are not",
          "intersecting in each fold)")

def report_label_distrib(labels):
    classes, counts = np.unique(labels, return_counts=True)
    norm_counts = counts/counts.sum()
    for i, label in enumerate(classes):
        print("class:", label, "norm count:", norm_counts[i])

def generate_partition_datatensor(pe_datatensor, data_partitions):
    datatensor_partitions = {}
    for run_num in data_partitions:
        datatensor_partitions[run_num] = {}
        for dsettype in data_partitions[run_num]:
            target_ids = data_partitions[run_num][dsettype]
            datatensor_partition = PartitionDataTensor(pe_datatensor, target_ids, dsettype, run_num)
            datatensor_partitions[run_num][dsettype] = datatensor_partition
    return(datatensor_partitions)

def construct_load_dataloaders(dataset_fold, dsettypes, config, wrk_dir):
    """construct dataloaders for the dataset for one fold

       Args:
            dataset_fold: dictionary,
                          example: {'train': <neural.dataset.PartitionDataTensor at 0x1cec95c96a0>,
                                    'validation': <neural.dataset.PartitionDataTensor at 0x1cec95c9208>,
                                    'test': <neural.dataset.PartitionDataTensor at 0x1cec95c9240>,
                                   }
            dsettype: list, ['train', 'validation', 'test']
            config: dict, {'batch_size': int, 'num_workers': int}
            wrk_dir: string, folder path
    """

    # setup data loaders
    data_loaders = {}
    epoch_loss_avgbatch = {}
    flog_out = {}
    score_dict = {}
    for dsettype in dsettypes:
        if(dsettype == 'train'):
            shuffle = True
            drop_last = True
        else:
            shuffle = False
            drop_last = False
        data_loaders[dsettype] = DataLoader(dataset_fold[dsettype],
                                            batch_size=config['batch_size'],
                                            shuffle=shuffle,
                                            num_workers=config['num_workers'],
                                            drop_last=drop_last)

        epoch_loss_avgbatch[dsettype] = []
        score_dict[dsettype] = ContModelScore(0, 0.0, 0.0) # best_epoch_indx, spearman_corr, pearson_corr 
        if wrk_dir:
            flog_out[dsettype] = os.path.join(wrk_dir, dsettype + ".log")
        else:
            flog_out[dsettype] = None

    return (data_loaders, epoch_loss_avgbatch, score_dict, flog_out)


class ConcatDataLoaders():
    """concatenate multiple dataloaders
    
    Args:
        dataloaders: list of torch.utils.DataLoaders
        mode: string, {'cycle', 'common_size'}. 
                      common_size: the loader will keep creating batches until it exhausts the smaller dataset and then exits
                      cycle: the loader will keep creating batches (cycles on smaller ones) until it exhausts the larger dataset and then exits
    
    
    """

    def __init__(self, dataloaders, dataset_names, mode='cycle'):
        self.dataloaders = dataloaders
        self.datasetnames = dataset_names
        self.mode = mode
        self.__len__()
        
    def __iter__(self):
        self.dloader_iterators = []
        for data_loader in self.dataloaders:
            if self.mode == 'cycle':
                if len(data_loader) < self.max_num_batches:
                    iterator = cycle(data_loader)
                else:
                    iterator = iter(data_loader)
            elif self.mode == 'common_size':
                iterator = iter(data_loader)
            self.dloader_iterators.append(iterator)
        return self

    def __next__(self):
        batches_lst = []
        for dloader_iter in self.dloader_iterators:
            batch = next(dloader_iter)
            batches_lst.append(batch)
        return tuple(batches_lst)

    def __len__(self):
        lst = sorted([len(dloader) for dloader in self.dataloaders])
        min_val = lst[0]
        max_val = lst[-1]
        # print('length of dloaders:', lst)
        # print('min:', min_val)
        # print('max:', max_val)
        # print()
        self.min_num_batches = min_val
        self.max_num_batches = max_val
        
        if self.mode == 'cycle':
            return max_val
        elif self.mode == 'common_size':
            return min_val

    
def construct_load_multiple_dataloaders(dataset_fold_lst, dsettypes, config, wrk_dir):
    """construct dataloaders for the dataset for one fold

       Args:
            dataset_fold_lst: list of dictionaries of the form,
                          example: {'train': <PartitionDataTensor at 0x1cec95c96a0>,
                                    'validation': <PartitionDataTensor at 0x1cec95c9208>,
                                    'test': <PartitionDataTensor at 0x1cec95c9240>,
                                   }
            dsettype: list, ['train', 'validation', 'test']
            config: dict, {'batch_size': int, 'num_workers': int}
            wrk_dir: string, folder path
    """

    # setup data loaders
    data_loaders = {}
    epoch_loss_avgbatch = {}
    flog_out = {}
    score_dict = {}
    print(config)
    loader_mode = config.get('loader_mode')
    datasets_name  = config.get('datasets_name')
    print('loader_mode:', loader_mode)
    print('datasets_name:', datasets_name)
    for dsettype in dsettypes:
        if(dsettype == 'train'):
            shuffle = True
            drop_last = True
        else:
            shuffle = False
            drop_last = False

        dloader_lst = []
        for tdict in dataset_fold_lst:
            dloader = DataLoader(tdict[dsettype],
                                batch_size=config['batch_size'],
                                shuffle=shuffle,
                                num_workers=config['num_workers'],
                                drop_last=drop_last)
            dloader_lst.append(dloader)
        data_loaders[dsettype] = ConcatDataLoaders(dloader_lst, datasets_name, mode=loader_mode)

        epoch_loss_avgbatch[dsettype] = []
        score_dict[dsettype] = ContModelScore(0, 0.0, 0.0) # best_epoch_indx, spearman_corr, pearson_corr 
        if wrk_dir:
            flog_out[dsettype] = os.path.join(wrk_dir, dsettype + ".log")
        else:
            flog_out[dsettype] = None

    return (data_loaders, epoch_loss_avgbatch, score_dict, flog_out)

def compute_weight(y):
    # weighting samples during loss computation as proposed by https://www.sciencedirect.com/science/article/pii/S0092867423003318
    
    bsize = y.shape[0]
    fill_val = -1.
    r = np.full((bsize, 2), fill_val)
    num = (np.array(y.cpu())+1.)**6
    denom = np.exp(18.)
    r[:,0] = 1. + num/denom
    r[:,1] = np.full((bsize,), 5.)
    return r.min(axis=1)

def compute_weight_23klib(y):    
    bsize = y.shape[0]
    fill_val = -1.
    r = np.full((bsize, 2), fill_val)
    num = (np.array(y.cpu())+1.)**6
    denom = np.exp(13.)
    r[:,0] = 1. + num/denom
    r[:,1] = np.full((bsize,), 5.)
    return r.min(axis=1)

def compute_correction_type_weight(weight_arr, correction_type_mat):
    #  'Correction_Deletion', 1, 3.81
    #  'Correction_Insertion', 2, 3.62
    #  'Correction_Replacement', 3, 2.17
    # [1, 0.7, 0.6]
    # np.array([[3.81, 3.62, 2.17]])
    return (weight_arr * correction_type_mat).sum(axis=1)