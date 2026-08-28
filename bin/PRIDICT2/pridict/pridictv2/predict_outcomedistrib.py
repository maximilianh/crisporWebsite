import os
import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.utils.data import DataLoader
from matplotlib import pyplot as plt
from tqdm import tqdm

from .hyperparam import get_saved_config
from .utilities import ReaderWriter, build_predictions_df, check_na,switch_layer_to_traineval_mode,require_nonleaf_grad
from .data_preprocess import PESeqProcessor
from .dataset import create_datatensor, MinMaxNormalizer
from ..rnn.rnn import RNN_Net

from .model import AnnotEmbeder_InitSeq, AnnotEmbeder_MutSeq, FeatureEmbAttention, \
                   MLPEmbedder, MLPDecoder, MaskGenerator, MLPDecoderDistribution
from .hyperparam import RNNHyperparamConfig
from .dataset import ConcatDataLoaders
from .data_preprocess import Viz_PESeqs, check_editing_alignment_correctness

class PRIEML_Model:
    def __init__(self, device, wsize=20, normalize='none', include_MFE=False, include_addendumfeat=False,  normalizer_dict=None, fdtype=torch.float32):
        self.device = device
        self.wsize = wsize
        self.normalize = normalize
        self.include_MFE = include_MFE
        self.include_addendumfeat = include_addendumfeat
        self.normalizer_dict = normalizer_dict
        self.fdtype = fdtype
        """
        base_90k: trained model using PRIDICT 1 schwank library
        base_390K: Multitask model trained on PRIDICT 1 schwank library and hyongbum ~300k library
        base_23k: multitask model trained on 23k library directly
        base_90k_decinit_HEKschwank_FT': fine tuned base_90k model on 23k library using HEK schwank decoder as initialization
        base_390k_decinit_HEKhyongbum_FT: fine tuned base_390k model on 23k library using HEK hyongbum decoder as initialization
        base_390k_decinit_HEKschwank_FT: fine tuned base_390k model on 23k library using HEK schwank decoder as initialization
        
        """
        self.modelnames_celltype_map = {'base_90k':['HEK'],
                                        'base_390k':['HEKschwank', 'HEKhyongbum'],
                                        'base_23k':['HEK', 'K562'],
                                        'base_90k_decinit_HEKschwank_FT':['HEK', 'K562'],
                                        'base_390k_decinit_HEKhyongbum_FT':['HEK', 'K562'],
                                        'base_390k_decinit_HEKschwank_FT':['HEK', 'K562']}
        
    def get_celltypes(self, modelname):
        return self.modelnames_celltype_map[modelname]

    def _process_df(self, df):
        """
        Args:
            df: pandas dataframe
        """
        print('--- processing input data frame ---')
        normalize = self.normalize
        assert normalize in {'none', 'max', 'minmax', 'standardize'}

        pe_seq_processor=PESeqProcessor()
        df = df.copy()
        # reset index in order to perform proper operations down the stream !
        df.reset_index(inplace=True, drop=True)
        if 'correction_type_categ' not in df:
            print('--- creating correction type categories ---')
            correction_categs = ['Deletion', 
                                'Insertion', 
                                'Replacement']
            df['correction_type_categ'] = pd.Categorical(df['Correction_Type'], categories=correction_categs)
            correction_type_df = pd.get_dummies(df['correction_type_categ'], prefix='Correction', prefix_sep='_')
            df = pd.concat([df, correction_type_df], axis=1)
        if 'seq_id' not in df:
            print('--- creating seq_id ---')
            df['seq_id'] = [f'seq_{i}' for i in range(df.shape[0])]

        # retrieve continuous column names
        minmax_normalizer = MinMaxNormalizer(include_MFE = self.include_MFE, 
                                             include_addendumfeat = self.include_addendumfeat,
                                             normalizer_dict=self.normalizer_dict)
        norm_colnames = minmax_normalizer.get_colnames()

        # print(norm_colnames)

        tdf, proc_seq_init_df,num_init_cols, proc_seq_mut_df, num_mut_cols = pe_seq_processor.process_init_mut_seqs(df,
                                                                                                               'wide_initial_target',
                                                                                                               'wide_mutated_target',
                                                                                                               align_symbol=2)
        df = pd.merge(left = df,
                      right = tdf[['seq_id', 
                                   'wide_initial_target_align', 
                                   'wide_mutated_target_align', 
                                   'Correction_Length_effective']],
                      how='inner',
                      left_on=['seq_id'],
                      right_on=['seq_id'])
        # print(df.columns.tolist())
        # add PBSlength in case it is part of colnames
        if 'PBSlength' in norm_colnames and 'PBSlength' not in df:
            print('--- creating PBSlength ---')
            df['PBSlength'] = proc_seq_init_df['end_PBS'] - proc_seq_init_df['start_PBS']
        if normalize != 'none':
            print('--- normalizing continuous features ---')
            norm_colnames = minmax_normalizer.normalize_cont_cols(df, normalize_opt=normalize, suffix='_norm')
        # print(norm_colnames)
        # make sure everything checks out
        check_na(proc_seq_init_df)
        check_na(proc_seq_mut_df)
        check_na(tdf)
        # print(check_editing_alignment_correctness(tdf, correction_len_colname='Correction_Length_effective'))
        return norm_colnames, df, proc_seq_init_df,num_init_cols, proc_seq_mut_df, num_mut_cols

    def _construct_datatensor(self, norm_colnames, df, proc_seq_init_df,num_init_cols, proc_seq_mut_df, num_mut_cols, y_ref=[]):
        print('--- creating datatensor ---')
        wsize=self.wsize # to read this from options dumped on disk
        dtensor = create_datatensor(df, 
                                    proc_seq_init_df, num_init_cols, 
                                    proc_seq_mut_df, num_mut_cols,
                                    norm_colnames,
                                    window=wsize, 
                                    y_ref=y_ref)
  
        return dtensor

    def _construct_dloader(self, dtensor, cell_types, batch_size):
        print('--- creating datatloader ---')
        dloader_lst = []
        for __ in cell_types:
            dloader = DataLoader(dtensor,
                                batch_size=batch_size,
                                shuffle=False,
                                drop_last=False,
                                num_workers=0,
                                sampler=None)
            dloader_lst.append(dloader)
        return ConcatDataLoaders(dloader_lst, cell_types, mode = 'cycle')

    def _load_model_config(self, mconfig_dir):
        print('--- loading model config ---')
        mconfig, options = get_saved_config(mconfig_dir)
        return mconfig, options

    def _build_base_model(self, config, cell_types=[]):
        print('--- building model ---')
        device = self.device
        mconfig, options = config
        model_config = mconfig['model_config']


        # print('mconfig:', mconfig)
        # print('options:', options)
        # print()
        
        if cell_types:
            datasets_name_lst = cell_types
        else:
            datasets_name_lst = options.get('datasets_name') # trained celltypes/datasets
            # fixes cases where name is composite and separated by _
            datasets_name_lst = ["".join(dname.split('_')) for dname in datasets_name_lst]
        # print('dataset_names_lst:', datasets_name_lst)


        separate_attention_layers = options.get('separate_attention_layers')
        separate_seqlevel_embedder = options.get('separate_seqlevel_embedder')

        # print('separate_seqlevel_embedder:', separate_seqlevel_embedder)


        fdtype = self.fdtype
        annot_embed = options.get('annot_embed')
        assemb_opt = options.get('assemb_opt')
        seqlevel_featdim = options.get('seqlevel_featdim')
        num_outcomes = options.get('num_outcomes', 3)
        loss_func = options.get('loss_func')

        init_annot_embed = AnnotEmbeder_InitSeq(embed_dim=model_config.embed_dim,
                                                annot_embed=annot_embed,
                                                assemb_opt=assemb_opt)
        mut_annot_embed = AnnotEmbeder_MutSeq(embed_dim=model_config.embed_dim,
                                              annot_embed=annot_embed,
                                              assemb_opt=assemb_opt)
        if assemb_opt == 'stack':
            init_embed_dim = model_config.embed_dim + 3*annot_embed
            mut_embed_dim = model_config.embed_dim + 2*annot_embed
            z_dim = np.min([init_embed_dim, mut_embed_dim])//2
        else:
            init_embed_dim = model_config.embed_dim
            mut_embed_dim = model_config.embed_dim
            z_dim = np.min([init_embed_dim, mut_embed_dim])//2 

        init_encoder = RNN_Net(input_dim =init_embed_dim,
                              hidden_dim=model_config.embed_dim,
                              z_dim=z_dim,
                              device=device,
                              num_hiddenlayers=model_config.num_hidden_layers,
                              bidirection=model_config.bidirection,
                              rnn_pdropout=model_config.p_dropout,
                              rnn_class=model_config.rnn_class,
                              nonlinear_func=model_config.nonlin_func,
                              fdtype=fdtype)
        mut_encoder= RNN_Net(input_dim =mut_embed_dim,
                              hidden_dim=model_config.embed_dim,
                              z_dim=z_dim,
                              device=device,
                              num_hiddenlayers=model_config.num_hidden_layers,
                              bidirection=model_config.bidirection,
                              rnn_pdropout=model_config.p_dropout,
                              rnn_class=model_config.rnn_class,
                              nonlinear_func=model_config.nonlin_func,
                              fdtype=fdtype)

        # seqlevel_featembeder = MLPEmbedder(inp_dim=seqlevel_featdim,
        #                                    embed_dim=z_dim,
        #                                    mlp_embed_factor=1,
        #                                    nonlin_func=model_config.nonlin_func,
        #                                    pdropout=model_config.p_dropout, 
        #                                    num_encoder_units=1)
        
        seqlevel_embedder_lst = []
        mlp_embed_factor = 1
        if separate_seqlevel_embedder:
            for dname in datasets_name_lst:
                seqlevel_featembeder = MLPEmbedder(inp_dim=seqlevel_featdim,
                                                embed_dim=z_dim,
                                                mlp_embed_factor=mlp_embed_factor,
                                                nonlin_func=model_config.nonlin_func,
                                                pdropout=model_config.p_dropout, 
                                                num_encoder_units=1)
                seqlevel_embedder_lst.append((seqlevel_featembeder, f'seqlevel_featembeder_{dname}'))
        else:
            seqlevel_featembeder = MLPEmbedder(inp_dim=seqlevel_featdim,
                                            embed_dim=z_dim,
                                            mlp_embed_factor=mlp_embed_factor,
                                            nonlin_func=model_config.nonlin_func,
                                            pdropout=model_config.p_dropout, 
                                            num_encoder_units=1)
            seqlevel_embedder_lst.append((seqlevel_featembeder, f'seqlevel_featembeder'))

        if separate_attention_layers:
            attn_modules_lst = []
            seq_types = ['init', 'mut']
            for dname in datasets_name_lst:
                tmp_lst = []
                for seq_type in seq_types: # original and mutated
                    for attn_type in ['local', 'global']:
                        tmp_lst.append((FeatureEmbAttention(z_dim), f'{attn_type}_featemb_{seq_type}_attn_{dname}'))
                attn_modules_lst.append(tmp_lst)
            # print('using separate attention layers!!')
            # print(attn_modules_lst)
        else:
            attn_modules_lst = []
            seq_types = ['init', 'mut']
            tmp_lst = []
            for seq_type in seq_types: # original and mutated
                for attn_type in ['local', 'global']:
                    tmp_lst.append((FeatureEmbAttention(z_dim), f'{attn_type}_featemb_{seq_type}_attn'))
            attn_modules_lst.append(tmp_lst)
       
        # print('loss_func:', loss_func)
        if loss_func in {'KLDloss', 'CEloss'}:
            # print('using MLPDecoderDistribution')
            decoder_lst = []
            for dname in datasets_name_lst:
                decoder  = MLPDecoderDistribution(5*z_dim,
                                                embed_dim=z_dim,
                                                outp_dim=num_outcomes,
                                                mlp_embed_factor=2,
                                                nonlin_func=model_config.nonlin_func, 
                                                pdropout=model_config.p_dropout, 
                                                num_encoder_units=1)
                decoder_lst.append((decoder,f'decoder_{dname}'))
        else:
            # print('using MLPDecoder')
            decoder_lst = []
            for dname in datasets_name_lst:
                decoder  = MLPDecoder(5*z_dim,
                                    embed_dim=z_dim,
                                    outp_dim=num_outcomes,
                                    mlp_embed_factor=2,
                                    nonlin_func=model_config.nonlin_func, 
                                    pdropout=model_config.p_dropout, 
                                    num_encoder_units=1,
                                    infer_sigma=False) 
                decoder_lst.append((decoder,f'decoder_{dname}'))


        # group submodels
        models = [(init_annot_embed, 'init_annot_embed'), 
                  (mut_annot_embed, 'mut_annot_embed'),
                  (init_encoder, 'init_encoder'),
                  (mut_encoder, 'mut_encoder')]
        
        for i_data in range(len(datasets_name_lst)):
            models += attn_modules_lst[i_data]

        models += seqlevel_embedder_lst
        models += decoder_lst
        
        # 4 + num_dsets*4 + num_dsest + num_dset
        return models

    def _load_model_statedict_(self, models, model_dir):
        print('--- loading trained model ---')
        device = self.device
        fdtype = self.fdtype
        # load state_dict pth
        state_dict_dir = os.path.join(model_dir, 'model_statedict')
        for m, m_name in models:
            m.load_state_dict(torch.load(os.path.join(state_dict_dir, '{}.pkl'.format(m_name)), map_location=device))
        # update model's fdtype and move to device
        for m, m_name in models:
            m.type(fdtype).to(device)
            m.eval()
        return models

    def _run_prediction(self, models, dloader, y_ref=[]):

        device = self.device
        fdtype = self.fdtype
        mask_gen = MaskGenerator()
        requires_grad = False
        pred_score = []

        num_loaders = len(dloader.dataloaders)
        
        #TODO: it assumes that all dataset in concatdataloader either have reference or not
        # make this flexible or check in different way 
        if dloader.dataloaders[0].dataset.y_score is not None:
            ref_score = []
        else:
            ref_score = None
        
        assert (len(y_ref) > 0 and len(y_ref) <= 3), '# of target outcomes should be > 0 and not exceed 3!'


        seqs_ids_lst = []

        # models = ((init_annot_embed, 'init_annot_embed'), 
        #           (mut_annot_embed, 'mut_annot_embed'),
        #           (init_encoder, 'init_encoder'),
        #           (mut_encoder, 'mut_encoder'),
        #           (local_featemb_init_attn, 'local_featemb_init_attn'),
        #           (local_featemb_mut_attn, 'local_featemb_mut_attn'),
        #           (global_featemb_init_attn, 'global_featemb_init_attn'),
        #           (global_featemb_mut_attn, 'global_featemb_mut_attn'),
        #           (seqlevel_featembeder, 'seqlevel_featembeder'),
        #           (decoder, 'decoder'))

        init_annot_embed = models[0][0]
        mut_annot_embed = models[1][0]
        init_encoder = models[2][0]
        mut_encoder = models[3][0]
    

        # local_featemb_init_attn = models[4][0]
        # local_featemb_mut_attn = models[5][0]
        # global_featemb_init_attn = models[6][0]
        # global_featemb_mut_attn = models[7][0]
        # decoder = models[9][0]
        datasets_name_lst = dloader.datasetnames
        offset_indx = 4
        # for i in range(num_loaders):
        #     dec_indx = offset_indx+(num_loaders*4)+(i+1)
        #     dec_name = models[dec_indx][-1].split('_')[-1]
        #     datasets_name_lst.append(dec_name)

        # print('datasets_name_lst:',datasets_name_lst)
        # for i, (m, mname) in enumerate(models):
        #     print('indx:', i, 'mname:', mname)
        #     print(m)
        #     print()
        

        if len(models) < 4 + len(datasets_name_lst)*4 + len(datasets_name_lst) + len(datasets_name_lst):
            one_seqlevelembedder_flag = True
        else:
            one_seqlevelembedder_flag = False


        dataset_ids_lst = []
        for i_batch, samples_batch_lst in tqdm(enumerate(dloader)):
            for i_data, samples_batch in enumerate(samples_batch_lst):
                X_init_nucl, X_init_proto, X_init_pbs, X_init_rt, \
                X_mut_nucl, X_mut_pbs, X_mut_rt, \
                x_init_len, x_mut_len, seqlevel_feat, \
                y_val, b_seqs_indx, b_seqs_id = samples_batch


                X_init_nucl = X_init_nucl.to(device)
                X_init_proto = X_init_proto.to(device)
                X_init_pbs = X_init_pbs.to(device)
                X_init_rt = X_init_rt.to(device)

                X_mut_nucl = X_mut_nucl.to(device)
                X_mut_pbs = X_mut_pbs.to(device)
                X_mut_rt = X_mut_rt.to(device)
                seqlevel_feat = seqlevel_feat.type(fdtype).to(device)
                
                
                with torch.set_grad_enabled(False):
                    X_init_batch = init_annot_embed(X_init_nucl, X_init_proto, X_init_pbs, X_init_rt)
                    X_mut_batch = mut_annot_embed(X_mut_nucl, X_mut_pbs, X_mut_rt)
                    # print('X_init_batch.shape:', X_init_batch.shape)
                    # print('X_mut_batch.shape:',X_mut_batch.shape)
                    # print('x_init_len.shape', x_init_len.shape)
                    # print('x_mut_len.shape:', x_mut_len.shape)

                    # print(np.unique(x_init_len))
                    # print(np.unique(x_mut_len))
                    # (bsize,)

                    
                    # masks (bsize, init_seqlen) or (bsize, mut_seqlen)
                    x_init_m = mask_gen.create_content_mask((X_init_batch.shape[0], X_init_batch.shape[1]), x_init_len)
                    x_mut_m =  mask_gen.create_content_mask((X_mut_batch.shape[0], X_mut_batch.shape[1]), x_mut_len)

                    __, z_init = init_encoder.forward_complete(X_init_batch, x_init_len, requires_grad=requires_grad)
                    __, z_mut =  mut_encoder.forward_complete(X_mut_batch, x_mut_len, requires_grad=requires_grad)
        
                    max_seg_len = z_init.shape[1]
                    init_mask = x_init_m[:,:max_seg_len].to(device)

                    # global attention
                    # s (bsize, embed_dim)

                    # print('accessing attention vectors:', offset_indx + (i_data*offset_indx)  + 0)
                    # print('model access:', offset_indx + (i_data*offset_indx)  + 0)
                    local_featemb_init_attn, __ = models[offset_indx + (i_data*offset_indx)  + 0]
                    # print(local_featemb_init_attn)
                    global_featemb_init_attn, __ = models[offset_indx + (i_data*offset_indx) + 1]
                    # print(global_featemb_init_attn)
                    local_featemb_mut_attn, __ = models[offset_indx + (i_data*offset_indx)   + 2]
                    # print(local_featemb_mut_attn)
                    global_featemb_mut_attn, __ = models[offset_indx + (i_data*offset_indx)  + 3]
                    # print(global_featemb_mut_attn)
    

                    s_init_global, __ = global_featemb_init_attn(z_init, mask=init_mask)
                    # local attention
                    s_init_local, __ = local_featemb_init_attn(z_init, mask=X_init_rt[:,:max_seg_len])

                    max_seg_len = z_mut.shape[1]
                    mut_mask = x_mut_m[:,:max_seg_len].to(device)
                    s_mut_global, __ = global_featemb_mut_attn(z_mut, mask=mut_mask)
                    s_mut_local, __ = local_featemb_mut_attn(z_mut, mask=X_mut_rt[:,:max_seg_len])
                     
                    # print('accessing seqlevelfeat_embedder:', offset_indx + (num_loaders*4) + i_data)
                    if one_seqlevelembedder_flag:
                        seqlevel_model_indx = 0
                        decoder_offset = 1
                    else:
                        seqlevel_model_indx = i_data
                        decoder_offset = num_loaders
                    seqlevel_featembeder, __ = models[offset_indx + (num_loaders*4) + seqlevel_model_indx]
                    # print('model access:', offset_indx + (num_loaders*4) + seqlevel_model_indx)

                    # print(seqlevel_featembeder)
                    seqfeat = seqlevel_featembeder(seqlevel_feat)
                    # y (bsize, 1)
                    # y_hat_logit, y_sigma = decoder(torch.cat([s_init_global, s_init_local, s_mut_global, s_mut_local, seqfeat], axis=-1))
                    
                    
                    decoder, __ = models[offset_indx + (num_loaders*4) + (decoder_offset) + i_data]

                    mlp_dec_flag = isinstance(decoder, MLPDecoder)

                    y_hat_logit = decoder(torch.cat([s_init_global, s_init_local, s_mut_global, s_mut_local, seqfeat], axis=-1))

                    num_outcomes = decoder.W_mu.out_features
                    # print('num_outcomes:', num_outcomes)
                    if mlp_dec_flag:
                        if ref_score is not None:
                            y_batch = y_val.type(fdtype) * 100.
                            ref_score.extend(y_batch[:,:num_outcomes].tolist())

                        pred_score.extend((torch.exp(y_hat_logit) - 1).tolist())
                        # print('torch.exp()-1')
                        # pred_score.extend((torch.exp(y_hat_logit-1)).tolist())


                    else:
                        if ref_score is not None:
                            y_batch = y_val.type(fdtype)
                            ref_score.extend(y_batch[:,:num_outcomes].tolist())

                        pred_score.extend(torch.exp(y_hat_logit).tolist())


                    seqs_ids_lst.extend(b_seqs_id.tolist())
                    dataset_ids_lst.extend([datasets_name_lst[i_data]]*len(b_seqs_id))

 
        predictions_df = build_predictions_df(seqs_ids_lst, 
                                              ref_score, 
                                              pred_score, 
                                              y_ref,
                                              dset_names=dataset_ids_lst)
        return predictions_df

    def prepare_data(self, df, model_name, cell_types=[], y_ref=[], batch_size=500):
        """
        Args:
            df: pandas dataframe
            y_ref: list (optional), list of reference outcome names
            batch_size: int, number of samples to process per batch
        """
        norm_colnames, df, proc_seq_init_df,num_init_cols, proc_seq_mut_df, num_mut_cols =  self._process_df(df)
        dtensor = self._construct_datatensor(norm_colnames, 
                                             df, 
                                             proc_seq_init_df,
                                             num_init_cols, 
                                             proc_seq_mut_df, 
                                             num_mut_cols, 
                                             y_ref=y_ref)
        
        # :class:`ConcatDataLoader`
        # TODO: to define use case where user limit what prediction heads to be used based on specified cell_types
        if cell_types:
            available_cell_types = cell_types
        else:
            available_cell_types = self.get_celltypes(model_name)



        dloader = self._construct_dloader(dtensor, available_cell_types, batch_size)
        return dloader


    def build_retrieve_models(self, model_dir):

        mconfig_dir = os.path.join(model_dir, 'config')
        mconfig = self._load_model_config(mconfig_dir)
        # print(mconfig)
        
        # list of tuples (model, model_name)
        models = self._build_base_model(mconfig)

        models = self._load_model_statedict_(models, model_dir)


        return models


    def predict_from_dloader(self, dloader, model_dir, y_ref=[]):
        
        cell_types = dloader.datasetnames
        mconfig_dir = os.path.join(model_dir, 'config')
        mconfig = self._load_model_config(mconfig_dir)
        
        # list of tuples (model, model_name)
        models = self._build_base_model(mconfig, cell_types)

        models = self._load_model_statedict_(models, model_dir)
        # print(f'running prediction for prime editor | model_dir: {model_dir}')
        pred_df = self._run_prediction(models, dloader, y_ref=y_ref)

        return pred_df

    def predict_from_dloader_using_loaded_models(self, dloader, models, y_ref=[]):
        pred_df = self._run_prediction(models, dloader, y_ref=y_ref)
        return pred_df

    def compute_avg_predictions(self, df):
        grp_cols = ['seq_id']
        if 'dataset_name' in df:
            grp_cols += ['dataset_name']
        agg_df = df.groupby(by=grp_cols).mean()
        agg_df.reset_index(inplace=True)
        for colname in ('run_num', 'Unnamed: 0'):
            if colname in agg_df:
                del agg_df[colname]
        return agg_df

    def visualize_seqs(self, df, seqsids_lst):
        """
        Args:
            df: pandas dataframe
            seqids_lst: list of sequence ids to plot

        """

        out_tb = {}
        tseqids = set(seqsids_lst).intersection(set(df['seq_id'].unique()))

        for tseqid in tqdm(tseqids):
            out_tb[tseqid] = Viz_PESeqs().viz_align_initmut_seq(df,
                                                                tseqid, 
                                                                wsize=self.wsize,
                                                                return_type='html')

        return out_tb