import re
import itertools
import os
import pandas as pd
import numpy as np
from scipy import stats
from prettytable import PrettyTable
from tqdm import tqdm
from Bio import pairwise2
from Bio.pairwise2 import format_alignment 

def get_char(seq):
    """split string int sequence of chars returned in pandas.Series"""
    chars = list(seq)
    return pd.Series(chars)

def is_initial_eq_mutated(df):
    check = (df['wide_initial_target'] == df['wide_mutated_target']).sum()
    if check:
        flag = True
    else:
        flag = False
    print('% of wide_target_initial == wide_muated_target:', 100*check/df.shape[0])
    return flag

def is_contiguous(lst, max_gap_length=1):
    forder_diff = np.diff(lst)
    max_val = forder_diff.max()
    min_val = forder_diff.min()
    gap_info = np.unique(forder_diff, return_counts=True)
    print('list of unique gaps:', gap_info[0])
    print('list of unique gaps distribution:',gap_info[1])
    print('max_gap:', max_val, 'min_gap:',min_val)
    flag = True
    try:
        assert max_val == min_val == max_gap_length
    except Exception as e:
        flag = False
    finally:
        return flag

def align_wt_mut_seqs_manual(sdf):

    if isinstance(sdf, pd.DataFrame):
        sdf = sdf.squeeze()

    correction_type = sdf['Correction_Type']
    correction_len = sdf['Correction_Length']
    # in case of {Insertion, Deletion}, we have only one starting position for these types of edits
    # case of {Replacement}, we have list of edited positions
    editpos_lst = get_startend_pos(sdf['deepeditposition_lst'])
    ewindow_st = editpos_lst[0]

    if correction_type in {'Insertion', 'Deletion'}:
        sdf['Correction_Length_effective'] = correction_len

    elif correction_type == 'Replacement':
        # we have multiple editing positions for this type of edit
        sdf['Correction_Length_effective'] = len(editpos_lst)
    
    # df['edit_pos_range'] = df['deepeditposition_lst'].str.strip('[]').str.split(',')
    # df['edit_pos_len'] = pd.Series([len(elm) for elm in df['edit_pos_range']])

    w = sdf['wide_initial_target']
    m = sdf['wide_mutated_target']

    # we just align when Insertion or Deletion occurs to save compute time
    if correction_type  == 'Insertion':
        w_repr = w[:ewindow_st] + '-'*correction_len + w[ewindow_st:]
        if w[ewindow_st:] != m[ewindow_st+correction_len:]:
            m_repr = m + '-'*correction_len 
        else:
            m_repr = m
    elif correction_type == 'Deletion':
        m_repr = m[:ewindow_st] + '-'*correction_len + m[ewindow_st:]
        if w[ewindow_st+correction_len:] != m[ewindow_st:]:
            w_repr = w + '-'*correction_len
        else:
            w_repr = w
    else:
        w_repr = w
        m_repr = m
    sdf['wide_initial_target_align'] = w_repr
    sdf['wide_mutated_target_align'] = m_repr
    return sdf

def align_wt_mut_seqs(sdf):
    
    if isinstance(sdf, pd.DataFrame):
        sdf = sdf.squeeze()

    correction_type = sdf['Correction_Type']
    correction_len = sdf['Correction_Length']
    # in case of {Insertion, Deletion}, we have only one starting position for these types of edits
    # case of {Replacement}, we have list of edited positions
    editpos_lst = get_startend_pos(sdf['deepeditposition_lst'])
    ewindow_st = editpos_lst[0]
    
    if correction_type in {'Insertion', 'Deletion'}:
        match = 1
        mismatch = -100
        gap_open = -1.
        gap_extend = -0.1
        penalize_end_gaps = False
        sdf['Correction_Length_effective'] = correction_len

    elif correction_type == 'Replacement':
        # we have multiple editing positions for this type of edit
        match = 1
        mismatch = 0
        gap_open = -100
        gap_extend = -100    
        penalize_end_gaps = True
        # effective number of changed bases
        sdf['Correction_Length_effective'] = len(editpos_lst)

        
    w = sdf['wide_initial_target']
    m = sdf['wide_mutated_target']
    
    # we just align when Insertion or Deletion occurs to save compute time
    if correction_type in {'Insertion', 'Deletion'}:
        alignments = pairwise2.align.globalms(w[ewindow_st:], m[ewindow_st:], 
                                              match, mismatch, gap_open, gap_extend,
                                              one_alignment_only=True,
                                              penalize_end_gaps=penalize_end_gaps)


        a,b,c=format_alignment(*alignments[0]).split('\n')[:-2]
        w_repr = w[:ewindow_st]+a
        m_repr = m[:ewindow_st]+c
    else:
        w_repr = w
        m_repr = m

    sdf['wide_initial_target_align'] = w_repr
    sdf['wide_mutated_target_align'] = m_repr

    return sdf

class PESeqProcessor:
    def __init__(self):
        """"""
        # target columns used/needed while processing the dataframe

        self.target_cols = ['seq_id',
                             'wide_initial_target',
                             'wide_mutated_target',
                             'deepeditposition',
                             'deepeditposition_lst', 
                             'Correction_Type', 
                             'Correction_Length',
                             'protospacerlocation_only_initial', 
                             'PBSlocation',
                             'RT_initial_location',
                             'RT_mutated_location']


    def align_seqs(self, df, init_seq_colname, mut_seq_colname):
        # harmonize the representation of sequences to UPPERCASE
        for colname in [init_seq_colname, mut_seq_colname]:
            df[colname] = df[colname].str.upper()
        # align sequences
        tqdm.pandas(desc='aligning sequences')
        a_df = df.groupby(by=['seq_id']).progress_apply(align_wt_mut_seqs_manual)
        a_df.reset_index(inplace=True, drop=True)
        r_df = pd.merge(left = df,
                        right = a_df[['seq_id', 'wide_initial_target_align', 
                                      'wide_mutated_target_align', 'Correction_Length_effective']],
                        how='inner',
                        left_on=['seq_id'],
                        right_on=['seq_id'])
        return r_df
    
    
    def mark_annotations(self, df):

        col_name_map = {'protospacerlocation_only_initial': 'Protos',
                        'PBSlocation': 'PBS',
                        'RT_initial_location':'RT_init',
                        'RT_mutated_location':'RT_mut'}

        df = df.copy()
        num_seqs = df.shape[0]

        tcols = ['protospacerlocation_only_initial', 
                 'PBSlocation',
                 'RT_initial_location',
                 'RT_mutated_location']
        # initial Insertion
        # mutated Deletion

        for colname in tqdm(tcols, total=len(tcols), desc='marking annotations'):

            range_df = df[colname].str.strip('[]').str.split(',')
            start_pos = range_df.str[0].values.astype(int)
            stop_pos =  range_df.str[-1].values.astype(int)

            df[f'start_{col_name_map[colname]}'] = start_pos
            df[f'end_{col_name_map[colname]}'] = stop_pos

            flag = False
            if colname == 'RT_initial_location':
                cond = df['Correction_Type'] == 'Insertion'
                flag = True
            elif colname == 'RT_mutated_location':
                cond = df['Correction_Type'] == 'Deletion'
                flag = True

            if flag:
                # adjusting the end of RT after alignment
                df.loc[cond, f'end_{col_name_map[colname]}'] +=  df.loc[cond, 'Correction_Length']

        # adjust the protospacer when the end position overlaps with Insertion cases in wide_initial_target
        cond = (df['Correction_Type'] == 'Insertion') & (df['deepeditposition'] < df['end_Protos'])
        df.loc[cond, f'end_Protos'] = df.loc[cond, f'end_Protos'] + df.loc[cond, 'Correction_Length']

        # create the EWindow marking for - symbol
        range_df = df['deepeditposition_lst'].str.strip('[]').str.split(',')
        start_pos = range_df.str[0].values.astype(int)
        stop_pos =  range_df.str[-1].values.astype(int)
        for seq_type in ('init', 'mut'):
            df[f'start_EWindow_{seq_type}'] = start_pos
            df[f'end_EWindow_{seq_type}'] = stop_pos

            if seq_type == 'init':
                cond = df['Correction_Type'] == 'Insertion'
            elif seq_type == 'mut':
                cond = df['Correction_Type'] == 'Deletion'

            df.loc[cond,f'end_EWindow_{seq_type}'] += df.loc[cond, 'Correction_Length']
            df.loc[~cond, f'end_EWindow_{seq_type}'] = df.loc[~cond, f'start_EWindow_{seq_type}']

        df['start_EWindow_Protos'] = df['end_Protos']
        df['end_EWindow_Protos'] = df['end_Protos']

        cond = (df['Correction_Type'] == 'Insertion') & (df['deepeditposition'] < df['end_Protos'])
        df.loc[cond, 'start_EWindow_Protos'] = df.loc[cond, 'deepeditposition']
        df.loc[cond, 'end_EWindow_Protos'] = df.loc[cond, f'start_EWindow_Protos'] + df.loc[cond, 'Correction_Length']

        return df 
 

    def add_seq_annotations(self, data_df, num_cols, align_symbol):
        """

        Args:
            proc_seq_df: pandas.DataFrame returned from :func:`self.process_perbase_df`
            data_df: pandas.DataFrame (dataset dataframe)
            num_cols: int, number of columns from identified nucleotides
            seq_type: string, in {initial, mutated}
        """
        data_df = data_df.copy()
        num_seqs = data_df.shape[0]

        tcolnames = ['Protos', 'PBS', 'RT_init', 'RT_mut']

        for colname in tqdm(tcolnames, total=len(tcolnames), desc='creating annotation matrices'):

            start_pos = data_df[f'start_{colname}']
            stop_pos =  data_df[f'end_{colname}']
            arr = np.full((num_seqs, num_cols), 0)

            if colname in {'RT_init', 'RT_mut', 'Protos'}:
                if colname == 'Protos':
                    suffix = 'Protos'
                else:
                    suffix = colname.split('_')[-1] # init or mut
                start_ewindow = data_df[f'start_EWindow_{suffix}']
                stop_ewindow  = data_df[f'end_EWindow_{suffix}']

                for i in tqdm(range(num_seqs), total=num_seqs, desc=f'adding annotation for {colname}'):
                    r = start_pos[i]
                    c = stop_pos[i]
                    e = start_ewindow[i]
                    s = stop_ewindow[i]
                    arr[i, r:e] = 1
                    arr[i, e:s] = align_symbol
                    arr[i, s:c] = 1
            else:
                for i in tqdm(range(num_seqs), total=num_seqs, desc=f'adding annotation for {colname}'):
                    r = start_pos[i]
                    c = stop_pos[i]
                    arr[i, r:c] = 1
            arr_df = pd.DataFrame(arr)
            arr_df.columns = [f'{colname}{i}' for i in range(0, num_cols)]
            data_df = pd.concat([data_df, arr_df], axis=1)

    #     # adjust protospacer when Insertion occurs (i.e. it has to extend when we add - in wide_initial_target_align)
    #     cond = (data_df['Correction_Type'] == 'Insertion') & (data_df['deepeditposition'] < data_df['end_Protos'])
    #     for indx in tindices:
    #         st_p = data_df.loc[indx, 'end_Protos']
    #         correction_len = t3.loc[indx, 'Correction_Length']
    #         t3.loc[f'Protos{i}' for i in range(st_p, st_p+correction_len)] = align_symb
    #         left_over = t3.loc[indx, 'end_Protos'] - t3.loc[indx, 'deepeditposition'] 
    #         t3.loc[f'Protos{i}' for i in range(st_p+correction_len, st_p+correction_len+left_over)] = 1

        # in this setup end_RT_mut == end_RT_init
        data_df['start_seq'] = data_df['start_Protos']
        data_df['end_seq'] = data_df['end_RT_mut']

        return data_df
    
    def split_matrices(self, df, num_cols):
        tcolnames = []
        for colname in ['Protos', 'PBS', 'RT_init']:
            tcolnames += [f'{colname}{i}' for i in range(0, num_cols)]
        for colname in ['Protos', 'PBS', 'RT_init', 'seq']:
            for prefix in ['start', 'end']:
                tcolnames += [f'{prefix}_{colname}']
        proc_init_df = df[tcolnames]
            
        tcolnames = []
        for colname in ['PBS', 'RT_mut']:
            tcolnames += [f'{colname}{i}' for i in range(0, num_cols)]
        for colname in ['PBS', 'RT_mut', 'seq']:
            for prefix in ['start', 'end']:
                tcolnames += [f'{prefix}_{colname}']
        proc_mut_df = df[tcolnames]
        return proc_init_df, proc_mut_df

    def process_init_mut_seq_visual(self, df, init_seq_colname, mut_seq_colname, align_symbol=2):
        """
        
        Args:
            df: pandas.DataFrame (dataset dataframe)
            init_seq_colname: string, column name of initial target sequence
            mut_seq_colname: string, column name of mutated target sequence

        """
        tdf = df[self.target_cols].copy()
        # this adds wide_initial_target_align and wide_mutated_target_align columns
        tdf = self.align_seqs(tdf, init_seq_colname, mut_seq_colname)
        # adds marking start and end of Protos, PBS, RT_init, RT_mut, and EWindows
        tdf = self.mark_annotations(tdf)
        # get maximum length of sequence
        # note since we align both sequences, we should have equal max length
        max_num_cols =  tdf[f'{init_seq_colname}_align'].str.len().max()
        # adds annotation columns Protos0 .. ProtosMaxNum, PBS0 .. PBSMaxNum, RT_init0 .. RT_initMaxNum, RT_mut0 .. RT_mutMaxNum,
        tdf = self.add_seq_annotations(tdf, max_num_cols, align_symbol)
        return tdf

    def process_init_mut_seqs(self, df, init_seq_colname, mut_seq_colname, align_symbol=2):
        """
        
        Args:
            df: pandas.DataFrame (dataset dataframe)
            init_seq_colname: string, column name of initial target sequence
            mut_seq_colname: string, column name of mutated target sequence

        """
        seqid_colname = 'seq_id'
        tdf = df[self.target_cols].copy()
        # this adds wide_initial_target_align and wide_mutated_target_align columns
        tdf = self.align_seqs(tdf, init_seq_colname, mut_seq_colname)
        # adds marking start and end of Protos, PBS, RT_init, RT_mut, and EWindows
        tdf = self.mark_annotations(tdf)
        # get maximum length of sequence
        # note since we align both sequences, we should have equal max length
        max_num_cols =  tdf[f'{init_seq_colname}_align'].str.len().max()
        # adds annotation columns Protos0 .. ProtosMaxNum, PBS0 .. PBSMaxNum, RT_init0 .. RT_initMaxNum, RT_mut0 .. RT_mutMaxNum,
        tdf = self.add_seq_annotations(tdf, max_num_cols, align_symbol)
        proc_init_df, proc_mut_df = self.split_matrices(tdf, max_num_cols)

        tqdm._instances.clear()
        pbar = tqdm(total=4)
        # adds the letter and encoding of the letter columns i.e. B0 .. BMaxNum, L0 .. LMaxNum, 
        init_df, num_init_cols = self.process_perbase_df(tdf, seqid_colname, f'{init_seq_colname}_align')
        pbar.update(1)
        mut_df, num_mut_cols = self.process_perbase_df(tdf, seqid_colname, f'{mut_seq_colname}_align')
        pbar.update(2)
        proc_seq_init_df = pd.concat([proc_init_df, init_df], axis=1)
        pbar.update(3)
        proc_seq_mut_df = pd.concat([proc_mut_df, mut_df], axis=1)
        pbar.update(4)
        pbar.close()
        self.validate_df(proc_seq_init_df)
        self.validate_df(proc_seq_mut_df)  

        return tdf, proc_seq_init_df, num_init_cols,  proc_seq_mut_df, num_mut_cols
    
        
    def process_perbase_df(self, df, seqid_colname, seq_colname):
        """cleans a data frame representing sequences and their edit info obtained from crispr experiment
        
        Args:
            df: pandas.DataFrame (dataset dataframe)
            target_cols: list of column names we need to keep
            seq_colname: string, column name of the target sequence (wild-type or mutated)
                    
        """
        pbar = tqdm(total=4)

        baseseq_df = df[seq_colname].apply(get_char)
        num_cols = baseseq_df.shape[1]
        baseseq_df.columns = [f'B{i}' for  i in range(0, num_cols)]
        pbar.update(2)

        baseseq_letters_df = baseseq_df.copy()
        baseseq_letters_df.columns = [f'L{i}' for  i in range(0, num_cols)]
        # replace base letters with numbers
        baseseq_df.replace(['A', 'C', 'T', 'G', '-', 'N'], [0,1,2,3,4,5], inplace=True)
        
        # replace Na in case of unequal length sequences
        baseseq_df = baseseq_df.fillna(5)
        baseseq_letters_df = baseseq_letters_df.fillna('N')
        pbar.update(3)

        base_df = pd.concat([df[[seqid_colname, seq_colname]],
                            baseseq_letters_df,
                            baseseq_df], axis=1)
        base_df.reset_index(inplace=True, drop=True)
        pbar.update(4)
        pbar.close()
        return base_df, num_cols
    
      
    def validate_df(self, df):
        assert df.isna().any().sum() == 0

def get_startend_pos(str_rpr):
    """parse string of format [pos1, pos2, etc..]"""
    lst = []
    for elm in str_rpr.strip('[]').split(','):
        lst.append(int(elm))
    return lst


class Viz_PESeqs:

    html_colors = {'blue':' #aed6f1',
                   'red':' #f5b7b1',
                   'green':' #a3e4d7',
                   'yellow':' #f9e79f',
                   'violet':'#d7bde2',
                   'brown_sugar':'#B47955',
                   'tan':'#DBB68F',
                   'tangerine':'#F19455',
                   'grey':'#BEBEBE'}
    
    codes = {'A':'@', 'C':'!', 'T':'#', 'G':'_', 'init':'~', 'mut':'%', 'ewindow':'§', 'prob':'`', '-':'/'}
    
    nucl_colrmap = {'A':'red',
                    'C':'yellow',
                    'T':'blue',
                    'G':'green',
                    '-':'grey',
                    'init':'violet',
                    'mut':'brown_sugar',
                    'ewindow':'tan',
                    'prob':'tangerine'}
    lower_thr = 0
    upper_thr = 99
    
    def __init__(self):
        pass

    @classmethod
    def viz_align_initmut_seq(clss, data_df, seqid, wsize=0, return_type='html'):
        """
        Args:
            data_df: dataset df 
            seqid: string, sequence id
            wsize: int, number of characters to include beyond the end of max([RT_initial, RT_mutated])
            return_type: string, default `html`
        """
        seqproc = PESeqProcessor()
        seqid_colname = 'seq_id'
        tdf = data_df.loc[data_df[seqid_colname] == seqid].copy()
        tdf = seqproc.process_init_mut_seq_visual(tdf, 'wide_initial_target', 'wide_mutated_target', align_symbol=2)
        return clss.viz_align_initmut_seq_precomputed(tdf, seqid, wsize=wsize, return_type=return_type)
        
    @classmethod
    def viz_align_initmut_seq_precomputed(clss, data_df, seqid, wsize=0, return_type='html'):
        """
        Args:
            data_df: dataset df 
            seqid: string, sequence id
            wsize: int, number of characters to include beyond the end of max([RT_initial, RT_mutated])
            return_type: string, default `html`
        """
        
        lower_thr = clss.lower_thr
        upper_thr = clss.upper_thr
        codes = clss.codes

        seqid_colname = 'seq_id'
        # squeeze turns pd.DataFrame. into pd.Series
        c_seq_datadf = data_df.loc[data_df[seqid_colname] == seqid].squeeze() 
        max_num_cols = len(c_seq_datadf['wide_initial_target_align'])
        print('max_num_cols:', max_num_cols)
        
        st_rt_mut, end_rt_mut = c_seq_datadf['start_RT_mut'],  c_seq_datadf['end_RT_mut']
        st_rt_wt, end_rt_wt = c_seq_datadf['start_RT_init'],  c_seq_datadf['end_RT_init']
        st_pbs, end_pbs = c_seq_datadf['start_PBS'],  c_seq_datadf['end_PBS']
        st_proto, end_proto =  c_seq_datadf['start_Protos'],  c_seq_datadf['end_Protos']
        annot_dict = {'Protos':(st_proto, end_proto),
                      'PBS':(st_pbs, end_pbs),
                      'RT_init':(st_rt_wt, end_rt_wt),
                      'RT_mut':(st_rt_mut, end_rt_mut)}
        print('annot_dict:',annot_dict)
        
        end_annot_seq = max([end_rt_mut, end_rt_wt])
        st_annot_seq =  st_proto
        print('end_annot_seq:', end_annot_seq)
        print('st_annot_seq:', st_annot_seq)


        max_upper_thr = np.min([st_annot_seq+upper_thr, end_annot_seq + wsize, max_num_cols])
        st_seq = np.clip(st_annot_seq - wsize, a_min=lower_thr, a_max=max_upper_thr)
        end_seq = np.clip(end_annot_seq + wsize, a_min=lower_thr, a_max=max_upper_thr)
        print('st_seq:', st_seq)
        print('end_seq:', end_seq)
        
        
        # align both sequences
        w_repr, m_repr = c_seq_datadf['wide_initial_target_align'],  c_seq_datadf['wide_mutated_target_align']

        
        tb = PrettyTable()
        tb.field_names = [f'Seq. ID:\n{seqid}'] + [f'{i}' for i in range(st_seq, end_seq)]
        
        # initial sequence information
        # Protos
        # PBS
        # RT
        # wide_initial_target
        if 'y' in c_seq_datadf:
            score = c_seq_datadf['y']
        else:
            score = None

        correction_type = c_seq_datadf['Correction_Type']
        correction_len = int(c_seq_datadf['Correction_Length'])
        editpos_lst = get_startend_pos(c_seq_datadf['deepeditposition_lst'])

        ewindow_st = editpos_lst[0]
        ewindow_end = end_seq
        print('ewindow_st:',ewindow_st)
        print('ewindow_end:', ewindow_end)

        
        if correction_type in {'Insertion', 'Deletion'}:
            ewindow_str_lst = [f'Editing: {correction_type}']
            tmp = ['' for elm in range(st_seq, ewindow_end)]
            offset = ewindow_st
            for i in range(correction_len):
                tmp[i+offset] = f"{codes['ewindow']}*"
            ewindow_str_lst += tmp
  
        elif correction_type == 'Replacement':
            ewindow_str_lst = [f'Editing: {correction_type}']
            tmp = ['' for elm in range(st_seq, ewindow_end)]
            for pos in editpos_lst:
                tmp[pos] = f"{codes['ewindow']}*"
            ewindow_str_lst += tmp
                
        # add annotations
        for colname in ('Protos', 'PBS', 'RT_init'):
            vals = c_seq_datadf[[f'{colname}{i}' for i in range(st_seq, end_seq)]].values
            print(vals)
            annot_lst = [f'{colname}']
            cl_lst = []
            for annot in vals:
                if annot == 0:
                    cl_lst += [' ']
                else:
                    cl_lst += [f"{codes['init']}{annot}"]
            annot_lst += cl_lst
            tb.add_row(annot_lst)

        if correction_type == 'Deletion':
            tb.add_row(ewindow_str_lst)
            
        wide_init_target = w_repr
        # print(w_repr, len(w_repr))
        # print(m_repr, len(m_repr))
        # print(st_seq, end_seq)
        init_target_str_lst = ['Initial sequence'] + [f'{codes[nucl]}{nucl}' for nucl in wide_init_target[st_seq:end_seq]]
        # print(init_target_str_lst, len(init_target_str_lst))
        tb.add_row(init_target_str_lst)
        
        # mutated sequence information
        # wide_mutated_target
        # PBS
        # RT
        
        end_mut_seq = end_seq
        wide_mut_target = m_repr
        if score is not None:
            mut_target_str_lst = ['{}Mutated sequence\n Prob. edit={:.3f}'.format(codes['prob'], score)] + \
                                [f'{codes[nucl]}{nucl}' for nucl in wide_mut_target[st_seq:end_seq]]
        else:
            mut_target_str_lst = ['{}Mutated sequence\n'.format(codes['prob'])] + \
                                 [f'{codes[nucl]}{nucl}' for nucl in wide_mut_target[st_seq:end_seq]] 
        
        tb.add_row(mut_target_str_lst)
        
        if correction_type in {'Insertion', 'Replacement'}:
            tb.add_row(ewindow_str_lst)
            
        for colname in ('PBS', 'RT_mut'):

            vals = c_seq_datadf[[f'{colname}{i}' for i in range(st_seq, end_seq)]].values

            print(vals)
            annot_lst = [f'{colname}']
            cl_lst = []
            for annot in vals:

                if annot == 0:
                    cl_lst += [' ']
                else:
                    cl_lst += [f"{codes['mut']}{annot}"]
            annot_lst += cl_lst
            tb.add_row(annot_lst)


        if return_type == 'html':
            return clss._format_html_table(tb.get_html_string())
        else: # default string
            return tb.get_string()


    @classmethod
    def _format_html_table(clss, html_str):
        html_colors = clss.html_colors
        codes = clss.codes
        nucl_colrmap = clss.nucl_colrmap
        for nucl in codes:
            ctext = codes[nucl]
            color = html_colors[nucl_colrmap[nucl]]
            html_str = re.sub(f'<td>{ctext}', '<td bgcolor="{}">'.format(color), html_str)
        return html_str

class PEVizFile:
    def __init__(self, resc_pth):
        # resc_pth: viz resources folder path
        # it contains 'header.txt',  'jupcellstyle.css', 'begin.txt', and 'end.txt'
        self.resc_pth = resc_pth
    def create(self, tablehtml, dest_pth, fname):
        resc_pth = self.resc_pth
        ls = []
        for ftname in ('header.txt',  'jupcellstyle.css', 'begin.txt'):
            with open(os.path.join(resc_pth, ftname), mode='r') as f:
                ls.extend(f.readlines())
        ls.append(tablehtml)
        with open(os.path.join(resc_pth, 'end.txt'), mode='r') as f:
            ls.extend(f.readlines())
        content = "".join(ls)
        with open(os.path.join(dest_pth, f'{fname}.html'), mode='w') as f:
            f.write(content)

def validate_df(df):
    print('number of NA:', df.isna().any().sum())

def check_editing_alignment_correctness(df, correction_len_colname='Correction_Length'):
    viol_seqs = []
    for gr_id, gr in tqdm(df.groupby(by=['seq_id'])):
        sdf = gr.squeeze()
        correction_type = sdf['Correction_Type']
        correction_len =  sdf[correction_len_colname]
        indx_positions = get_startend_pos(sdf['deepeditposition_lst'])
        indx_pos = indx_positions[0]
        
            
        w_repr = sdf['wide_initial_target_align']
        m_repr = sdf['wide_mutated_target_align']
        
        try:
            # assert that the sequences are equivalent before reaching the editing position
            note = 'edit_position'
            assert w_repr[:indx_pos] == m_repr[:indx_pos]

            if correction_type == 'Insertion':
                # XXXXX-----XXXXXXXXX------
                # XXXXXIIIIIXXXXXXXXX------
                padding_str = w_repr[indx_pos:indx_pos+correction_len]
                note = 'correction_len'
                assert set(padding_str) == {'-'}

                assert w_repr[indx_pos+correction_len:len(w_repr)-correction_len] == m_repr[indx_pos+correction_len:len(m_repr)-correction_len]

            elif correction_type == 'Deletion':
                # XXXXXIIIIIXXXXXXXXX------
                # XXXXX-----XXXXXXXXX------
                padding_str = m_repr[indx_pos:indx_pos+correction_len]
                note = 'correction_len'
                assert set(padding_str) == {'-'}

                assert m_repr[indx_pos+correction_len:len(m_repr)-correction_len] == w_repr[indx_pos+correction_len:len(w_repr)-correction_len]

            elif correction_type == 'Replacement':
                # XXXXXRRXXRXXXXXX
                # XXXXXXXXXXXXXXXX
                note = 'correction_len'
                assert len(indx_positions) == correction_len 
                pos_set = set(indx_positions)
                for i in range(len(w_repr)):
                    if i in pos_set:
                        note = 'edit_pos_!='
                        assert w_repr[i]!=m_repr[i]
                        
                    else:
                        note = 'edit_pos_=='
                        assert w_repr[i]==m_repr[i]
                        
            else:
                print('we should not be here')
        except Exception as e:
#             print(e)
#             print(gr_id, correction_type)
            viol_seqs.append((gr_id, correction_type, correction_len, note))
#             viol_seqs.append(gr)
#     viol_seqs_df = pd.concat(viol_seqs, ignore_index=True)
    if len(viol_seqs):
        viol_seqs_df = pd.DataFrame(viol_seqs)
        viol_seqs_df.columns = ['seq_id', 'correction_type', 'correction_len', 'issue']
    else:
        viol_seqs_df = []
    return viol_seqs_df

