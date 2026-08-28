import os
import itertools
import numpy as np
import torch
from torch import nn
import torch.multiprocessing as mp
from tqdm import tqdm
from .utilities import get_device, create_directory, ReaderWriter, dump_dict_content, \
                       perfmetric_report_cont, \
                       plot_loss, plot_xy, build_predictions_df,\
                       restrict_grad_, get_trainable_params,get_num_trainable_params,freeze_layers, \
                       perfmetric_report_multidata_cont
from .model import AnnotEmbeder_InitSeq, AnnotEmbeder_MutSeq, SH_Attention, FeatureEmbAttention, \
                   MLPEmbedder, MLPDecoder, MaskGenerator, MLPDecoderDistribution, init_params_
from ..rnn.rnn import RNN_Net

from .dataset import construct_load_dataloaders, construct_load_multiple_dataloaders, compute_correction_type_weight
from .hyperparam import RNNHyperparamConfig,  get_hyperparam_options
from .loss import CELoss, BalancedMSELoss



def generate_models_config(hyperparam_config, 
                           experiment_options, 
                           loss_func='RMSEloss'):
    
    dataloader_config = {'batch_size': hyperparam_config.batch_size,
                         'num_workers': 0}

    config = {'dataloader_config': dataloader_config,
              'model_config': hyperparam_config,
             }

    # experiment_options = {'experiment_desc',
    #                       'model_name', 
    #                       'input_dim',
    #                       'run_num', 
    #                       'max_timedelta_encoding', 
    #                       'fdtype', 
    #                       'y_col'}

    options = {'num_epochs': hyperparam_config.num_epochs,
               'weight_decay': hyperparam_config.l2_reg,
               'to_gpu':True,
               'loss_func':loss_func
               }

    options.update(experiment_options)

    return config, options

def build_config_map(trf_tup, experiment_options, loss_func='NPloss'):
    model_name = experiment_options['model_name']

    if(model_name == 'PE_RNN'):
        hyperparam_config = RNNHyperparamConfig(*trf_tup)
    elif(model_name == 'PE_RNN_distribution'):
        hyperparam_config = RNNHyperparamConfig(*trf_tup)
    elif(model_name == 'PE_RNN_distribution_multidata'):
        hyperparam_config = RNNHyperparamConfig(*trf_tup)
    run_num = -1 
    fdtype = torch.float32
    experiment_options['run_num'] = run_num
    experiment_options['fdtype'] = fdtype
    mconfig, options = generate_models_config(hyperparam_config,
                                              experiment_options,
                                              loss_func=loss_func)
    return mconfig, options

# def run_cont_pe_RNN(data_partition, dsettypes, config, options, wrk_dir, state_dict_dir=None, to_gpu=True, gpu_index=0):
#     pid = "{}".format(os.getpid())  # process id description
#     # get data loader config
#     dataloader_config = config['dataloader_config']
#     cld = construct_load_dataloaders(data_partition, dsettypes, dataloader_config, wrk_dir)
#     # dictionaries by dsettypes
#     data_loaders, epoch_loss_avgbatch, score_dict, flog_out = cld

#     device = get_device(to_gpu, gpu_index)  # gpu device
#     fdtype = options['fdtype']
#     loss_type = options.get('loss_func', 'mse')

#     if loss_type == 'MSEloss' or loss_type == 'RMSEloss':
#         loss_func = nn.MSELoss(reduction='mean')
#     elif loss_type == 'L1loss':
#         loss_func = nn.L1Loss(reduction='mean')
#     elif loss_type == 'Huberloss':
#         loss_func = nn.SmoothL1Loss(reduction='mean')

#     num_epochs = options.get('num_epochs', 50)
#     run_num = options.get('run_num')
#     # parse config dict
#     model_config = config['model_config']
#     model_name = options['model_name']
    
#     target_names = ['averageedited', 'averageunedited', 'averageindel']
#     target_outcome = options.get('target_outcome')
#     outcome_indx = target_names.index(target_outcome)

#     if(model_name == 'PE_RNN'):
#         annot_embed = options.get('annot_embed')
#         assemb_opt = options.get('assemb_opt')
#         seqlevel_featdim = options.get('seqlevel_featdim')

#         init_annot_embed = AnnotEmbeder_InitSeq(embed_dim=model_config.embed_dim,
#                                                 annot_embed=annot_embed,
#                                                 assemb_opt=assemb_opt)
#         mut_annot_embed = AnnotEmbeder_MutSeq(embed_dim=model_config.embed_dim,
#                                               annot_embed=annot_embed,
#                                               assemb_opt=assemb_opt)
#         if assemb_opt == 'stack':
#             init_embed_dim = model_config.embed_dim + 3*annot_embed
#             mut_embed_dim = model_config.embed_dim + 2*annot_embed
#             z_dim = np.min([init_embed_dim, mut_embed_dim])//2
#         else:
#             init_embed_dim = model_config.embed_dim
#             mut_embed_dim = model_config.embed_dim
#             z_dim = np.min([init_embed_dim, mut_embed_dim])//2 

#         init_encoder = RNN_Net(input_dim =init_embed_dim,
#                               hidden_dim=model_config.embed_dim,
#                               z_dim=z_dim,
#                               device=device,
#                               num_hiddenlayers=model_config.num_hidden_layers,
#                               bidirection=model_config.bidirection,
#                               rnn_pdropout=model_config.p_dropout,
#                               rnn_class=model_config.rnn_class,
#                               nonlinear_func=model_config.nonlin_func,
#                               fdtype=fdtype)
#         mut_encoder= RNN_Net(input_dim =mut_embed_dim,
#                               hidden_dim=model_config.embed_dim,
#                               z_dim=z_dim,
#                               device=device,
#                               num_hiddenlayers=model_config.num_hidden_layers,
#                               bidirection=model_config.bidirection,
#                               rnn_pdropout=model_config.p_dropout,
#                               rnn_class=model_config.rnn_class,
#                               nonlinear_func=model_config.nonlin_func,
#                               fdtype=fdtype)

#         local_featemb_init_attn = FeatureEmbAttention(z_dim)
#         local_featemb_mut_attn = FeatureEmbAttention(z_dim)

#         global_featemb_init_attn = FeatureEmbAttention(z_dim)
#         global_featemb_mut_attn = FeatureEmbAttention(z_dim)

#         seqlevel_featembeder = MLPEmbedder(inp_dim=seqlevel_featdim,
#                                            embed_dim=z_dim,
#                                            mlp_embed_factor=1,
#                                            nonlin_func=model_config.nonlin_func,
#                                            pdropout=model_config.p_dropout, 
#                                            num_encoder_units=1)

#         decoder  = MLPDecoder(5*z_dim,
#                               embed_dim=z_dim,
#                               outp_dim=1,
#                               mlp_embed_factor=2,
#                               nonlin_func=model_config.nonlin_func, 
#                               pdropout=model_config.p_dropout, 
#                               num_encoder_units=1)

#         # define optimizer and group parameters
#         models = ((init_annot_embed, 'init_annot_embed'), 
#                   (mut_annot_embed, 'mut_annot_embed'),
#                   (init_encoder, 'init_encoder'),
#                   (mut_encoder, 'mut_encoder'),
#                   (local_featemb_init_attn, 'local_featemb_init_attn'),
#                   (local_featemb_mut_attn, 'local_featemb_mut_attn'),
#                   (global_featemb_init_attn, 'global_featemb_init_attn'),
#                   (global_featemb_mut_attn, 'global_featemb_mut_attn'),
#                   (seqlevel_featembeder, 'seqlevel_featembeder'),
#                   (decoder, 'decoder'))
#         models_param = []
#         for m, m_name in models:
#             models_param.extend(list(m.parameters()))
#             init_params_(m)

#     if(state_dict_dir):  # load state dictionary of saved models
#         for m, m_name in models:
#             m.load_state_dict(torch.load(os.path.join(state_dict_dir, '{}.pkl'.format(m_name)), map_location=device))

#     # update models fdtype and move to device
#     for m, m_name in models:
#         m.type(fdtype).to(device)

#     if('train' in data_loaders):
#         weight_decay = options.get('weight_decay', 1e-4)
#         # see paper Cyclical Learning rates for Training Neural Networks for parameters' choice
#         # `https://arxive.org/pdf/1506.01186.pdf`
#         # pytorch version >1.1, scheduler should be called after optimizer
#         # for cyclical lr scheduler, it should be called after each batch update
#         num_iter = len(data_loaders['train'])  # num_train_samples/batch_size
#         c_step_size = int(np.ceil(5*num_iter))  # this should be 2-10 times num_iter
#         base_lr = 3e-4
#         max_lr = 5*base_lr  # 3-5 times base_lr
#         optimizer = torch.optim.Adam(models_param, weight_decay=weight_decay, lr=base_lr)
#         cyc_scheduler = torch.optim.lr_scheduler.CyclicLR(optimizer, 
#                                                           base_lr, 
#                                                           max_lr, 
#                                                           step_size_up=c_step_size,
#                                                           mode='triangular', 
#                                                           cycle_momentum=False)


#     if ('validation' in data_loaders):
#         m_state_dict_dir = create_directory(os.path.join(wrk_dir, 'model_statedict'))

#     if(num_epochs > 1):
#         fig_dir = create_directory(os.path.join(wrk_dir, 'figures'))

#     # dump config dictionaries on disk
#     config_dir = create_directory(os.path.join(wrk_dir, 'config'))
#     ReaderWriter.dump_data(config, os.path.join(config_dir, 'mconfig.pkl'))
#     ReaderWriter.dump_data(options, os.path.join(config_dir, 'exp_options.pkl'))
#     # store attention weights for validation and test set
#     seqid_fattnw_map = {dsettype: {} for dsettype in data_loaders if dsettype in {'test'}}
#     mask_gen = MaskGenerator()
#     print('updated!!')
#     for epoch in range(num_epochs):
#         # print("-"*35)
#         for dsettype in dsettypes:
#             print("device: {} | experiment_desc: {} | run_num: {} | epoch: {} | dsettype: {} | pid: {}"
#                   "".format(device, options.get('experiment_desc'), run_num, epoch, dsettype, pid))
#             pred_score = []
#             ref_score = []
#             seqs_ids_lst = []
#             data_loader = data_loaders[dsettype]

#             epoch_loss = 0.

#             if(dsettype == 'train'):  # should be only for train
#                 for m, m_name in models:
#                     m.train()
#                 requires_grad=True
#             else:
#                 for m, m_name in models:
#                     m.eval()
#                 requires_grad = False

#             for i_batch, samples_batch in enumerate(tqdm(data_loader)):
#                 # print('batch num:', i_batch)

#                 X_init_nucl, X_init_proto, X_init_pbs, X_init_rt, \
#                 X_mut_nucl, X_mut_pbs, X_mut_rt, \
#                 x_init_len, x_mut_len, seqlevel_feat, \
#                 y_val, b_seqs_indx, b_seqs_id = samples_batch


#                 X_init_nucl = X_init_nucl.to(device)
#                 X_init_proto = X_init_proto.to(device)
#                 X_init_pbs = X_init_pbs.to(device)
#                 X_init_rt = X_init_rt.to(device)

#                 X_mut_nucl = X_mut_nucl.to(device)
#                 X_mut_pbs = X_mut_pbs.to(device)
#                 X_mut_rt = X_mut_rt.to(device)
#                 seqlevel_feat = seqlevel_feat.type(fdtype).to(device)
                
#                 with torch.set_grad_enabled(dsettype == 'train'):
#                     X_init_batch = init_annot_embed(X_init_nucl, X_init_proto, X_init_pbs, X_init_rt)
#                     X_mut_batch = mut_annot_embed(X_mut_nucl, X_mut_pbs, X_mut_rt)
#                     # print('X_init_batch.shape:', X_init_batch.shape)
#                     # print('X_mut_batch.shape:',X_mut_batch.shape)
#                     # print('x_init_len.shape', x_init_len.shape)
#                     # print('x_mut_len.shape:', x_mut_len.shape)

#                     # print(np.unique(x_init_len))
#                     # print(np.unique(x_mut_len))
#                     # (bsize,)
#                     y_batch = 100*y_val.type(fdtype).to(device)[:, outcome_indx]
                    
#                     # masks (bsize, init_seqlen) or (bsize, mut_seqlen)
#                     x_init_m = mask_gen.create_content_mask((X_init_batch.shape[0], X_init_batch.shape[1]), x_init_len)
#                     x_mut_m =  mask_gen.create_content_mask((X_mut_batch.shape[0], X_mut_batch.shape[1]), x_mut_len)

#                     __, z_init = init_encoder.forward_complete(X_init_batch, x_init_len, requires_grad=requires_grad)
#                     __, z_mut =  mut_encoder.forward_complete(X_mut_batch, x_mut_len, requires_grad=requires_grad)
        
#                     max_seg_len = z_init.shape[1]
#                     init_mask = x_init_m[:,:max_seg_len].to(device)
#                     # global attention
#                     # s (bsize, embed_dim)
#                     s_init_global, __ = global_featemb_init_attn(z_init, mask=init_mask)
#                     # local attention
#                     s_init_local, __ = local_featemb_init_attn(z_init, mask=X_init_rt[:,:max_seg_len])

#                     max_seg_len = z_mut.shape[1]
#                     mut_mask = x_mut_m[:,:max_seg_len].to(device)
#                     s_mut_global, __ = global_featemb_mut_attn(z_mut, mask=mut_mask)
#                     s_mut_local, __ = local_featemb_mut_attn(z_mut, mask=X_mut_rt[:,:max_seg_len])

#                     seqfeat = seqlevel_featembeder(seqlevel_feat)
#                     # y (bsize, 1)
#                     # y_hat_logit, y_sigma = decoder(torch.cat([s_init_global, s_init_local, s_mut_global, s_mut_local, seqfeat], axis=-1))
                    
#                     y_hat_logit = decoder(torch.cat([s_init_global, s_init_local, s_mut_global, s_mut_local, seqfeat], axis=-1))

#                     loss = loss_func(y_batch, y_hat_logit)

#                     # if loss_type != 'Logploss':
#                     #     loss = loss_func(y_batch, y_hat_logit)
#                     # else:
#                     #     y_dist = torch.distributions.normal.Normal(y_hat_logit, y_sigma)
#                     #     loss = loss_func(y_batch, y_dist)

#                     if(dsettype == 'train'):
#                         # zero gradient
#                         optimizer.zero_grad()
#                         # print("computing loss")
#                         # backward step (i.e. compute gradients)
#                         loss.backward()
#                         # optimzer step -- update weights
#                         optimizer.step()
#                         # after each batch step the scheduler
#                         cyc_scheduler.step()

#                         # print(optimizer)

#                     pred_score.extend(y_hat_logit.view(-1).tolist())
#                     ref_score.extend(y_batch.view(-1).tolist())
#                     seqs_ids_lst.extend(list(b_seqs_id))

#                     epoch_loss += loss.item()

#             # end of epoch
#             epoch_loss_avgbatch[dsettype].append(epoch_loss/len(data_loader))
#             modelscore = perfmetric_report_cont(pred_score, ref_score, 
#                                                 epoch_loss_avgbatch[dsettype][-1], 
#                                                 epoch, flog_out[dsettype])
#             perf = modelscore.spearman_corr
#             if(perf > score_dict[dsettype].spearman_corr):
#                 score_dict[dsettype] = modelscore
#                 if(dsettype == 'validation'):
#                     for m, m_name in models:
#                         torch.save(m.state_dict(), os.path.join(m_state_dict_dir, '{}.pkl'.format(m_name)))
#                 if dsettype in {'test', 'validation'}:
#                     predictions_df = build_predictions_df(seqs_ids_lst, ref_score, pred_score, target_names)
#                     predictions_path = os.path.join(wrk_dir, f'predictions_{dsettype}.csv')
#                     predictions_df.to_csv(predictions_path)

#     if(num_epochs > 1):
#         for dsettype in epoch_loss_avgbatch:
#             plot_xy(np.arange(num_epochs), 
#                     epoch_loss_avgbatch[dsettype], 
#                     'number of epochs', 
#                     f'{loss_type} loss', 
#                     'epoch batch average loss', 
#                     dsettype,
#                     fig_dir)
#     # dump_scores
#     dump_dict_content(score_dict, list(score_dict.keys()), 'score', wrk_dir)


def run_cont_pe_RNN_distribution(data_partition, dsettypes, config, options, wrk_dir, state_dict_dir=None, to_gpu=True, gpu_index=0):
    pid = "{}".format(os.getpid())  # process id description
    # get data loader config
    dataloader_config = config['dataloader_config']
    cld = construct_load_dataloaders(data_partition, dsettypes, dataloader_config, wrk_dir)
    # dictionaries by dsettypes
    data_loaders, epoch_loss_avgbatch, score_dict, flog_out = cld

    device = get_device(to_gpu, gpu_index)  # gpu device
    fdtype = options['fdtype']
    loss_type = options.get('loss_func')

    if loss_type == 'MSEloss' or loss_type == 'RMSEloss':
        loss_func = nn.MSELoss(reduction='mean')
    elif loss_type == 'L1loss':
        loss_func = nn.L1Loss(reduction='mean')
    elif loss_type == 'Huberloss':
        loss_func = nn.SmoothL1Loss(reduction='mean')
    elif loss_type == 'KLDloss':
        loss_func = nn.KLDivLoss(reduction='none')
    elif loss_type == 'CEloss':
        loss_func = CELoss(reduction='none')
    
    print('loss_type:', loss_type)

    num_epochs = options.get('num_epochs', 50)
    run_num = options.get('run_num')
    # parse config dict
    model_config = config['model_config']
    model_name = options['model_name']
    num_outcomes = options['num_outcomes']

    target_names = ['averageedited', 'averageunedited', 'averageindel']

     

    if(model_name == 'PE_RNN_distribution'):
        annot_embed = options.get('annot_embed')
        assemb_opt = options.get('assemb_opt')
        seqlevel_featdim = options.get('seqlevel_featdim')

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

        local_featemb_init_attn = FeatureEmbAttention(z_dim)
        local_featemb_mut_attn = FeatureEmbAttention(z_dim)

        global_featemb_init_attn = FeatureEmbAttention(z_dim)
        global_featemb_mut_attn = FeatureEmbAttention(z_dim)

        seqlevel_featembeder = MLPEmbedder(inp_dim=seqlevel_featdim,
                                           embed_dim=z_dim,
                                           mlp_embed_factor=1,
                                           nonlin_func=model_config.nonlin_func,
                                           pdropout=model_config.p_dropout, 
                                           num_encoder_units=1)

        decoder  = MLPDecoderDistribution(5*z_dim,
                                        embed_dim=z_dim,
                                        outp_dim=num_outcomes,
                                        mlp_embed_factor=2,
                                        nonlin_func=model_config.nonlin_func, 
                                        pdropout=model_config.p_dropout, 
                                        num_encoder_units=1)

        # define optimizer and group parameters
        models = ((init_annot_embed, 'init_annot_embed'), 
                  (mut_annot_embed, 'mut_annot_embed'),
                  (init_encoder, 'init_encoder'),
                  (mut_encoder, 'mut_encoder'),
                  (local_featemb_init_attn, 'local_featemb_init_attn'),
                  (local_featemb_mut_attn, 'local_featemb_mut_attn'),
                  (global_featemb_init_attn, 'global_featemb_init_attn'),
                  (global_featemb_mut_attn, 'global_featemb_mut_attn'),
                  (seqlevel_featembeder, 'seqlevel_featembeder'),
                  (decoder, 'decoder'))
        models_param = []
        for m, m_name in models:
            models_param.extend(list(m.parameters()))
            init_params_(m)

    if(state_dict_dir):  # load state dictionary of saved models
        for m, m_name in models:
            m.load_state_dict(torch.load(os.path.join(state_dict_dir, '{}.pkl'.format(m_name)), map_location=device))

    # update models fdtype and move to device
    for m, m_name in models:
        m.type(fdtype).to(device)

    if('train' in data_loaders):
        weight_decay = options.get('weight_decay', 1e-4)
        # see paper Cyclical Learning rates for Training Neural Networks for parameters' choice
        # `https://arxive.org/pdf/1506.01186.pdf`
        # pytorch version >1.1, scheduler should be called after optimizer
        # for cyclical lr scheduler, it should be called after each batch update
        num_iter = len(data_loaders['train'])  # num_train_samples/batch_size
        c_step_size = int(np.ceil(5*num_iter))  # this should be 2-10 times num_iter
        base_lr = 3e-4
        max_lr = 5*base_lr  # 3-5 times base_lr
        optimizer = torch.optim.Adam(models_param, weight_decay=weight_decay, lr=base_lr)
        cyc_scheduler = torch.optim.lr_scheduler.CyclicLR(optimizer, 
                                                          base_lr, 
                                                          max_lr, 
                                                          step_size_up=c_step_size,
                                                          mode='triangular', 
                                                          cycle_momentum=False)


    if ('validation' in data_loaders):
        m_state_dict_dir = create_directory(os.path.join(wrk_dir, 'model_statedict'))

    if(num_epochs > 1):
        fig_dir = create_directory(os.path.join(wrk_dir, 'figures'))

    # dump config dictionaries on disk
    config_dir = create_directory(os.path.join(wrk_dir, 'config'))
    ReaderWriter.dump_data(config, os.path.join(config_dir, 'mconfig.pkl'))
    ReaderWriter.dump_data(options, os.path.join(config_dir, 'exp_options.pkl'))
    # store attention weights for validation and test set
    seqid_fattnw_map = {dsettype: {} for dsettype in data_loaders if dsettype in {'test'}}
    mask_gen = MaskGenerator()

    for epoch in range(num_epochs):
        # print("-"*35)
        for dsettype in dsettypes:
            print("device: {} | experiment_desc: {} | run_num: {} | epoch: {} | dsettype: {} | pid: {}"
                  "".format(device, options.get('experiment_desc'), run_num, epoch, dsettype, pid))
            pred_score = []
            ref_score = []
            seqs_ids_lst = []
            data_loader = data_loaders[dsettype]

            epoch_loss = 0.

            if(dsettype == 'train'):  # should be only for train
                for m, m_name in models:
                    m.train()
                requires_grad=True
            else:
                for m, m_name in models:
                    m.eval()
                requires_grad = False

            for i_batch, samples_batch in enumerate(tqdm(data_loader)):
                # print('batch num:', i_batch)

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
                
                
                with torch.set_grad_enabled(dsettype == 'train'):
                    X_init_batch = init_annot_embed(X_init_nucl, X_init_proto, X_init_pbs, X_init_rt)
                    X_mut_batch = mut_annot_embed(X_mut_nucl, X_mut_pbs, X_mut_rt)
                    # print('X_init_batch.shape:', X_init_batch.shape)
                    # print('X_mut_batch.shape:',X_mut_batch.shape)
                    # print('x_init_len.shape', x_init_len.shape)
                    # print('x_mut_len.shape:', x_mut_len.shape)

                    # print(np.unique(x_init_len))
                    # print(np.unique(x_mut_len))
                    # (bsize,)
                    y_batch = y_val.type(fdtype).to(device)
                    
                    # masks (bsize, init_seqlen) or (bsize, mut_seqlen)
                    x_init_m = mask_gen.create_content_mask((X_init_batch.shape[0], X_init_batch.shape[1]), x_init_len)
                    x_mut_m =  mask_gen.create_content_mask((X_mut_batch.shape[0], X_mut_batch.shape[1]), x_mut_len)

                    __, z_init = init_encoder.forward_complete(X_init_batch, x_init_len, requires_grad=requires_grad)
                    __, z_mut =  mut_encoder.forward_complete(X_mut_batch, x_mut_len, requires_grad=requires_grad)
        
                    max_seg_len = z_init.shape[1]
                    init_mask = x_init_m[:,:max_seg_len].to(device)
                    # global attention
                    # s (bsize, embed_dim)
                    s_init_global, __ = global_featemb_init_attn(z_init, mask=init_mask)
                    # local attention
                    s_init_local, __ = local_featemb_init_attn(z_init, mask=X_init_rt[:,:max_seg_len])

                    max_seg_len = z_mut.shape[1]
                    mut_mask = x_mut_m[:,:max_seg_len].to(device)
                    s_mut_global, __ = global_featemb_mut_attn(z_mut, mask=mut_mask)
                    s_mut_local, __ = local_featemb_mut_attn(z_mut, mask=X_mut_rt[:,:max_seg_len])

                    seqfeat = seqlevel_featembeder(seqlevel_feat)
                    # y (bsize, 1)
                    # y_hat_logit, y_sigma = decoder(torch.cat([s_init_global, s_init_local, s_mut_global, s_mut_local, seqfeat], axis=-1))
                    
                    y_hat_logit = decoder(torch.cat([s_init_global, s_init_local, s_mut_global, s_mut_local, seqfeat], axis=-1))

                    loss = loss_func(y_hat_logit, y_batch)
                    # print('loss.shape:',loss.shape)
                    loss = loss.sum(axis=-1).mean()
                    # print('mean_loss:', loss)

                    # if loss_type != 'Logploss':
                    #     loss = loss_func(y_batch, y_hat_logit)
                    # else:
                    #     y_dist = torch.distributions.normal.Normal(y_hat_logit, y_sigma)
                    #     loss = loss_func(y_batch, y_dist)

                    if(dsettype == 'train'):
                        optimizer.zero_grad()
                        # print("computing loss")
                        # backward step (i.e. compute gradients)
                        loss.backward()
                        # optimzer step -- update weights
                        optimizer.step()
                        # after each batch step the scheduler
                        cyc_scheduler.step()

                        # print(optimizer)

                    pred_score.extend(torch.exp(y_hat_logit).tolist())
                    ref_score.extend(y_batch.tolist())
                    seqs_ids_lst.extend(list(b_seqs_id))

                    epoch_loss += loss.item()

            # end of epoch
            epoch_loss_avgbatch[dsettype].append(epoch_loss/len(data_loader))
            modelscore = perfmetric_report_cont(pred_score, ref_score, 
                                                epoch_loss_avgbatch[dsettype][-1], 
                                                epoch, flog_out[dsettype])
            perf = modelscore.spearman_corr
            if(perf > score_dict[dsettype].spearman_corr):
                score_dict[dsettype] = modelscore
                if(dsettype == 'validation'):
                    for m, m_name in models:
                        torch.save(m.state_dict(), os.path.join(m_state_dict_dir, '{}.pkl'.format(m_name)))
                    ReaderWriter.dump_data({'epoch':epoch+1}, os.path.join(m_state_dict_dir, 'best_epoch.pkl'))
                if dsettype in {'test', 'validation'}:
                    predictions_df = build_predictions_df(seqs_ids_lst, ref_score, pred_score, target_names)
                    predictions_path = os.path.join(wrk_dir, f'predictions_{dsettype}.csv')
                    predictions_df.to_csv(predictions_path)

    if(num_epochs > 1):
        for dsettype in epoch_loss_avgbatch:
            plot_xy(np.arange(num_epochs), 
                    epoch_loss_avgbatch[dsettype], 
                    'number of epochs', 
                    f'{loss_type} loss', 
                    'epoch batch average loss', 
                    dsettype,
                    fig_dir)
    # dump_scores
    dump_dict_content(score_dict, list(score_dict.keys()), 'score', wrk_dir)


def run_cont_pe_RNN_distribution_multidata(data_partition, dsettypes, config, options, wrk_dir, 
                                                 state_dict_dir=None, to_gpu=True, gpu_index=0):
    pid = "{}".format(os.getpid())  # process id description

    # random_seed = options['random_seed']
    # torch.manual_seed(random_seed)
    # torch.cuda.manual_seed(random_seed)
    # torch.cuda.manual_seed_all(random_seed)
    # np.random.seed(random_seed)

    # get data loader config
    dataloader_config = config['dataloader_config']
    dataloader_config['loader_mode'] = options.get('loader_mode')
    dataloader_config['datasets_name'] = options.get('datasets_name')
    cld = construct_load_multiple_dataloaders(data_partition, dsettypes, dataloader_config, wrk_dir)
    # dictionaries by dsettypes
    data_loaders, epoch_loss_avgbatch, score_dict, flog_out = cld

    device = get_device(to_gpu, gpu_index)  # gpu device
    fdtype = options['fdtype']
    loss_type = options.get('loss_func')

    if loss_type == 'MSEloss' or loss_type == 'RMSEloss':
        loss_func = nn.MSELoss(reduction='none')
    elif loss_type == 'L1loss':
        loss_func = nn.L1Loss(reduction='none')
    elif loss_type == 'Huberloss':
        loss_func = nn.SmoothL1Loss(reduction='none')
    elif loss_type == 'KLDloss':
        loss_func = nn.KLDivLoss(reduction='none')
    elif loss_type == 'CEloss':
        loss_func = CELoss(reduction='none')
    elif loss_type == 'BalancedMSEloss':
        loss_func = BalancedMSELoss()

    print('loss_type:', loss_type)

    num_epochs = options.get('num_epochs', 50)
    run_num = options.get('run_num')
    # parse config dict
    model_config = config['model_config']
    model_name = options['model_name']
    num_outcomes = options['num_outcomes']
        
    # default_outcomes = ['averageedited', 'averageunedited', 'averageindel']

    target_names = options.get('target_names')

    datasets_name_lst = options.get('datasets_name')
    separate_attention_layers = options.get('separate_attention_layers')
    separate_seqlevel_embedder = options.get('separate_seqlevel_embedder')
    
    print('datasets_name_lst:', datasets_name_lst)

    weight_func_pointers = options.get('weight_func_pointers', [None for __ in range(len(datasets_name_lst))])
    correctiontype_weights = options.get('correctiontype_weights', [None for __ in range(len(datasets_name_lst))])
    
    if(model_name == 'PE_RNN_distribution_multidata'):
        annot_embed = options.get('annot_embed')
        assemb_opt = options.get('assemb_opt')
        seqlevel_featdim = options.get('seqlevel_featdim')

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

        # local_featemb_init_attn = FeatureEmbAttention(z_dim)
        # local_featemb_mut_attn = FeatureEmbAttention(z_dim)

        # global_featemb_init_attn = FeatureEmbAttention(z_dim)
        # global_featemb_mut_attn = FeatureEmbAttention(z_dim)
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
            print('using separate attention layers!!')
            print(attn_modules_lst)
        else:
            attn_modules_lst = []
            seq_types = ['init', 'mut']
            tmp_lst = []
            for seq_type in seq_types: # original and mutated
                for attn_type in ['local', 'global']:
                    tmp_lst.append((FeatureEmbAttention(z_dim), f'{attn_type}_featemb_{seq_type}_attn'))
            attn_modules_lst.append(tmp_lst)

        if loss_type in {'KLDloss', 'CEloss'}:
            print('using MLPDecoderDistribution')
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
            print('using MLPDecoder')
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




        # define optimizer and group parameters
        models = [(init_annot_embed, 'init_annot_embed'), 
                  (mut_annot_embed, 'mut_annot_embed'),
                  (init_encoder, 'init_encoder'),
                  (mut_encoder, 'mut_encoder')]
        models += seqlevel_embedder_lst
        
        for i_data in range(len(datasets_name_lst)):
            models += attn_modules_lst[i_data]

        models += decoder_lst
        print(models)

        models_param = []
        for m, m_name in models:
            models_param.extend(list(m.parameters()))
            init_params_(m)

    if(state_dict_dir):  # load state dictionary of saved models
        for m, m_name in models:
            m.load_state_dict(torch.load(os.path.join(state_dict_dir, '{}.pkl'.format(m_name)), map_location=device))

    # update models fdtype and move to device
    for m, m_name in models:
        m.type(fdtype).to(device)

    if('train' in data_loaders):
        weight_decay = options.get('weight_decay', 1e-4)
        # see paper Cyclical Learning rates for Training Neural Networks for parameters' choice
        # `https://arxive.org/pdf/1506.01186.pdf`
        # pytorch version >1.1, scheduler should be called after optimizer
        # for cyclical lr scheduler, it should be called after each batch update
        num_iter = len(data_loaders['train'])  # num_train_samples/batch_size
        c_step_size = int(np.ceil(5*num_iter))  # this should be 2-10 times num_iter
        base_lr = 3e-4
        max_lr = 5*base_lr  # 3-5 times base_lr
        optimizer = torch.optim.Adam(models_param, weight_decay=weight_decay, lr=base_lr)
        cyc_scheduler = torch.optim.lr_scheduler.CyclicLR(optimizer, 
                                                          base_lr, 
                                                          max_lr, 
                                                          step_size_up=c_step_size,
                                                          mode='triangular', 
                                                          cycle_momentum=False)


    if ('validation' in data_loaders):
        m_state_dict_dir = create_directory(os.path.join(wrk_dir, 'model_statedict'))

    if(num_epochs > 1):
        fig_dir = create_directory(os.path.join(wrk_dir, 'figures'))

    # dump config dictionaries on disk
    config_dir = create_directory(os.path.join(wrk_dir, 'config'))
    ReaderWriter.dump_data(config, os.path.join(config_dir, 'mconfig.pkl'))
    ReaderWriter.dump_data(options, os.path.join(config_dir, 'exp_options.pkl'))
    # store attention weights for validation and test set
    seqid_fattnw_map = {dsettype: {} for dsettype in data_loaders if dsettype in {'test'}}
    mask_gen = MaskGenerator()

    for epoch in range(num_epochs):
        # print("-"*35)
        for dsettype in dsettypes:
            print("device: {} | experiment_desc: {} | run_num: {} | epoch: {} | dsettype: {} | pid: {}"
                  "".format(device, options.get('experiment_desc'), run_num, epoch, dsettype, pid))
            pred_score = []
            ref_score = []
            seqs_ids_lst = []
            dataset_ids_lst = []
            
            data_loader = data_loaders[dsettype]

            epoch_loss = 0.

            if(dsettype == 'train'):  # should be only for train
                for m, m_name in models:
                    m.train()
                requires_grad=True
            else:
                for m, m_name in models:
                    m.eval()
                requires_grad = False

            for i_batch, samples_batch_lst in enumerate(tqdm(data_loader)):
                loss_tensor = torch.zeros(len(samples_batch_lst)).type(fdtype).to(device)
                for i_data, samples_batch in enumerate(samples_batch_lst):
                    # print('batch num:', i_batch)

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
                    
                    
                    with torch.set_grad_enabled(dsettype == 'train'):
                        X_init_batch = init_annot_embed(X_init_nucl, X_init_proto, X_init_pbs, X_init_rt)
                        X_mut_batch = mut_annot_embed(X_mut_nucl, X_mut_pbs, X_mut_rt)
                        # print('X_init_batch.shape:', X_init_batch.shape)
                        # print('X_mut_batch.shape:',X_mut_batch.shape)
                        # print('x_init_len.shape', x_init_len.shape)
                        # print('x_mut_len.shape:', x_mut_len.shape)

                        # print(np.unique(x_init_len))
                        # print(np.unique(x_mut_len))
                        # (bsize,)
                        y_batch = y_val.type(fdtype).to(device)
                        
                        # masks (bsize, init_seqlen) or (bsize, mut_seqlen)
                        x_init_m = mask_gen.create_content_mask((X_init_batch.shape[0], X_init_batch.shape[1]), x_init_len)
                        x_mut_m =  mask_gen.create_content_mask((X_mut_batch.shape[0], X_mut_batch.shape[1]), x_mut_len)

                        __, z_init = init_encoder.forward_complete(X_init_batch, x_init_len, requires_grad=requires_grad)
                        __, z_mut =  mut_encoder.forward_complete(X_mut_batch, x_mut_len, requires_grad=requires_grad)
            
                        max_seg_len = z_init.shape[1]
                        init_mask = x_init_m[:,:max_seg_len].to(device)

                        # global attention
                        # s (bsize, embed_dim)
                        # assign each attention layer
                        if separate_attention_layers:
                            attn_module_indx = i_data
                        else:
                            attn_module_indx = 0

                        local_featemb_init_attn, __ = attn_modules_lst[attn_module_indx][0]
                        global_featemb_init_attn, __ = attn_modules_lst[attn_module_indx][1]
                        local_featemb_mut_attn, __ = attn_modules_lst[attn_module_indx][2]
                        global_featemb_mut_attn, __ = attn_modules_lst[attn_module_indx][3]
        

                        s_init_global, __ = global_featemb_init_attn(z_init, mask=init_mask)
                        # local attention
                        s_init_local, __ = local_featemb_init_attn(z_init, mask=X_init_rt[:,:max_seg_len])

                        max_seg_len = z_mut.shape[1]
                        mut_mask = x_mut_m[:,:max_seg_len].to(device)
                        s_mut_global, __ = global_featemb_mut_attn(z_mut, mask=mut_mask)
                        s_mut_local, __ = local_featemb_mut_attn(z_mut, mask=X_mut_rt[:,:max_seg_len])
                        
                        seqlevel_featembeder, __ = seqlevel_embedder_lst[i_data]
                        seqfeat = seqlevel_featembeder(seqlevel_feat)
                        # y (bsize, 1)
                        # y_hat_logit, y_sigma = decoder(torch.cat([s_init_global, s_init_local, s_mut_global, s_mut_local, seqfeat], axis=-1))
                        
                        decoder, __ = decoder_lst[i_data]
                        y_hat_logit = decoder(torch.cat([s_init_global, s_init_local, s_mut_global, s_mut_local, seqfeat], axis=-1))
                        

                        if loss_type  == 'BalancedMSEloss':
                            y_batch = y_batch * 100.
                            # transforms scores from 0-100 to 0-4.x
                            y_batch_transf = torch.log1p(y_batch[:,:num_outcomes])
                
                            pred_score.extend((torch.exp(y_hat_logit) - 1).tolist())
                            ref_score.extend(y_batch[:,:num_outcomes].tolist())
                            loss = loss_func(y_hat_logit, y_batch_transf,  seqlevel_feat[:,1:4])

                        elif loss_type not in {'KLDloss', 'CEloss', 'BalancedMSEloss'}:
                            y_batch = y_batch * 100.
                            # transforms scores from 0-100 to 0-4.x
                            y_batch_transf = torch.log1p(y_batch[:,:num_outcomes])
                            y_weights = y_batch[:, 0] # first column is averageediting
                
                            pred_score.extend((torch.exp(y_hat_logit) - 1).tolist())
                            ref_score.extend(y_batch[:,:num_outcomes].tolist())

                            loss = loss_func(y_hat_logit, y_batch_transf)

                        elif loss_type in {'KLDloss', 'CEloss'}: # case of KLDloss or CEloss
                            y_batch_transf = y_batch[:,:num_outcomes]
                            y_weights = y_batch[:, 0] * 100. # first column is averageediting

                            pred_score.extend(torch.exp(y_hat_logit).tolist())
                            ref_score.extend(y_batch[:,:num_outcomes].tolist())

                            loss = loss_func(y_hat_logit, y_batch_transf)


                        # print('y_hat_logit.shape:', y_hat_logit.shape)
                        # print('y_batch_transf.shape:', y_batch_transf.shape)
                        # print('loss.shape:',loss.shape)
                        if loss_type == 'BalancedMSEloss':
                            loss_tensor[i_data] = loss
                        else:
                            weight_func = weight_func_pointers[i_data]
                            if weight_func is not None:
                                weights = weight_func(y_weights)
                                weights = torch.from_numpy(weights).type(fdtype).to(device).reshape(-1,1)
                                # print('weights.shape:', weights.shape)
                                wloss = loss * weights
                                # print('wloss.shape:', wloss.shape)
                            else:
                                wloss = loss
                                # print('wloss.shape:', wloss.shape)

                        
                        # # amplify losses based on correction types
                        # ctp_weight = correctiontype_weights[i_data]
                        # # print('ctp_weight:', ctp_weight)
                        # if ctp_weight is not None:
                        #    # seqlevel_feat[:,1:4] are 'Correction_Deletion', 'Correction_Insertion', 'Correction_Replacement
                        #    correction_type_mat = seqlevel_feat[:,1:4].cpu().numpy()
                        #    val = compute_correction_type_weight(ctp_weight, correction_type_mat)
                        # #    print('ctp_weight.shape: ', val.shape)
                        #    wloss *= torch.from_numpy(val).type(fdtype).to(device).reshape(-1,1)
                        # #    print('wloss.shape:', wloss.shape)

                            if loss_type in {'KLDloss', 'CEloss'}:
                                wloss = wloss.sum(axis=-1)
                            # print('wloss.shape:', wloss.shape)
                            loss_tensor[i_data] = wloss.mean(axis=0) # average across batches
                            

                        # pred_score.extend(torch.exp(y_hat_logit).tolist())
                        # ref_score.extend(y_batch.tolist())
                        seqs_ids_lst.extend(list(b_seqs_id))
                        dataset_ids_lst.extend([datasets_name_lst[i_data]]*len(b_seqs_id))

                    # print('mean_loss:', loss)

                    # if loss_type != 'Logploss':
                    #     loss = loss_func(y_batch, y_hat_logit)
                    # else:
                    #     y_dist = torch.distributions.normal.Normal(y_hat_logit, y_sigma)
                    #     loss = loss_func(y_batch, y_dist)
                
                loss = loss_tensor.sum()

                if(dsettype == 'train'):
                    optimizer.zero_grad()
                    # print("computing loss")
                    # backward step (i.e. compute gradients)
                    # print('loss_tensor:', loss_tensor)
                    loss.backward()
                    # optimzer step -- update weights
                    optimizer.step()
                    # after each batch step the scheduler
                    cyc_scheduler.step()

                    # print(optimizer)

                    
                epoch_loss += loss.item()

            # end of epoch
            epoch_loss_avgbatch[dsettype].append(epoch_loss/len(data_loader))
            # get average modelscore across the datasets
            modelscore = perfmetric_report_multidata_cont(pred_score, 
                                                          ref_score, 
                                                          dataset_ids_lst,
                                                          seqs_ids_lst,
                                                          epoch_loss_avgbatch[dsettype][-1], 
                                                          epoch, flog_out[dsettype])
            
            perf = modelscore.spearman_corr
            # print('modelscore.spearman_corr:', modelscore.spearman_corr)
            # print('modelscore.pearson_corr:', modelscore.pearson_corr)
            # print('modelscore.spearman_lst:', modelscore.spearman_lst)
            # print('modelscore.pearson_lst:', modelscore.pearson_lst)
            # print('modelscore.modelscores_lst:', modelscore.modelscores_lst)

            if(perf > score_dict[dsettype].spearman_corr):
                # print(f'score_dict[{dsettype}].spearman_corr:', score_dict[dsettype].spearman_corr)

                score_dict[dsettype] = modelscore
                
                if(dsettype == 'validation'):
                    for m, m_name in models:
                        torch.save(m.state_dict(), os.path.join(m_state_dict_dir, '{}.pkl'.format(m_name)))
                    ReaderWriter.dump_data({'epoch':epoch+1}, os.path.join(m_state_dict_dir, 'best_epoch.pkl'))
                if dsettype in {'test', 'validation'}:
                    predictions_df = build_predictions_df(seqs_ids_lst, 
                                                          ref_score, 
                                                          pred_score, 
                                                          target_names,
                                                          dset_names=dataset_ids_lst)
                    predictions_path = os.path.join(wrk_dir, f'predictions_{dsettype}.csv')
                    predictions_df.to_csv(predictions_path)

            # if(epoch == 9 or dsettype == 'test'):
            #     # print(f'score_dict[{dsettype}].spearman_corr:', score_dict[dsettype].spearman_corr)

            #     score_dict[dsettype] = modelscore
                
            #     if(dsettype == 'validation'):
            #         for m, m_name in models:
            #             torch.save(m.state_dict(), os.path.join(m_state_dict_dir, '{}.pkl'.format(m_name)))
            #         ReaderWriter.dump_data({'epoch':epoch+1}, os.path.join(m_state_dict_dir, 'best_epoch.pkl'))
            #     if dsettype in {'test', 'validation'}:
            #         predictions_df = build_predictions_df(seqs_ids_lst, 
            #                                               ref_score, 
            #                                               pred_score, 
            #                                               target_names,
            #                                               dset_names=dataset_ids_lst)
            #         predictions_path = os.path.join(wrk_dir, f'predictions_{dsettype}.csv')
            #         predictions_df.to_csv(predictions_path)

    if(num_epochs > 1):
        for dsettype in epoch_loss_avgbatch:
            plot_xy(np.arange(num_epochs), 
                    epoch_loss_avgbatch[dsettype], 
                    'number of epochs', 
                    f'{loss_type} loss', 
                    'epoch batch average loss', 
                    dsettype,
                    fig_dir)
    # dump_scores
    dump_dict_content(score_dict, list(score_dict.keys()), 'score', wrk_dir)

def run_tune_cont_pe_RNN_distribution_multidata(data_partition, dsettypes, config, options, wrk_dir, 
                                                 state_dict_dir=None, to_gpu=True, gpu_index=0):
    pid = "{}".format(os.getpid())  # process id description
    # get data loader config
    dataloader_config = config['dataloader_config']
    dataloader_config['loader_mode'] = options.get('loader_mode')
    dataloader_config['datasets_name'] = options.get('datasets_name')
    cld = construct_load_multiple_dataloaders(data_partition, dsettypes, dataloader_config, wrk_dir)
    # dictionaries by dsettypes
    data_loaders, epoch_loss_avgbatch, score_dict, flog_out = cld

    device = get_device(to_gpu, gpu_index)  # gpu device
    fdtype = options['fdtype']
    loss_type = options.get('loss_func')
    base_model_suffix = options.get('base_model_suffix')

    if loss_type == 'MSEloss' or loss_type == 'RMSEloss':
        loss_func = nn.MSELoss(reduction='mean')
    elif loss_type == 'L1loss':
        loss_func = nn.L1Loss(reduction='mean')
    elif loss_type == 'Huberloss':
        loss_func = nn.SmoothL1Loss(reduction='mean')
    elif loss_type == 'KLDloss':
        loss_func = nn.KLDivLoss(reduction='none')
    elif loss_type == 'CEloss':
        loss_func = CELoss(reduction='none')
        print('loss_type:', loss_type)

    num_epochs = options.get('num_epochs', 50)
    run_num = options.get('run_num')
    # parse config dict
    model_config = config['model_config']
    model_name = options['model_name']
    num_outcomes = options['num_outcomes']

    target_names = options.get('target_names')
    datasets_name_lst = options.get('datasets_name')
    separate_attention_layers = options.get('separate_attention_layers')
    separate_seqlevel_embedder = options.get('separate_seqlevel_embedder')
    
    print('datasets_name_lst:', datasets_name_lst)


    trainable_layernames = options['trainable_layernames']

    if(model_name == 'PE_RNN_distribution_multidata'):
        annot_embed = options.get('annot_embed')
        assemb_opt = options.get('assemb_opt')
        seqlevel_featdim = options.get('seqlevel_featdim')

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

        # local_featemb_init_attn = FeatureEmbAttention(z_dim)
        # local_featemb_mut_attn = FeatureEmbAttention(z_dim)

        # global_featemb_init_attn = FeatureEmbAttention(z_dim)
        # global_featemb_mut_attn = FeatureEmbAttention(z_dim)

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
            print('using separate attention layers!!')
            print(attn_modules_lst)
        else:
            attn_modules_lst = []
            seq_types = ['init', 'mut']
            tmp_lst = []
            for seq_type in seq_types: # original and mutated
                for attn_type in ['local', 'global']:
                    tmp_lst.append((FeatureEmbAttention(z_dim), f'{attn_type}_featemb_{seq_type}_attn'))
            attn_modules_lst.append(tmp_lst)

        if loss_type in {'KLDloss', 'CEloss'}:
            print('using MLPDecoderDistribution')
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
            print('using MLPDecoder')
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

      # define optimizer and group parameters
        models = [(init_annot_embed, 'init_annot_embed'), 
                  (mut_annot_embed, 'mut_annot_embed'),
                  (init_encoder, 'init_encoder'),
                  (mut_encoder, 'mut_encoder')]
        models += seqlevel_embedder_lst
        
        for i_data in range(len(datasets_name_lst)):
            models += attn_modules_lst[i_data]

        models += decoder_lst
        print(models)


        models_param = []
        for m, m_name in models:
            if m_name in trainable_layernames:
                models_param.extend(list(m.parameters()))

    if(state_dict_dir):  # load state dictionary of saved models
        for m, m_name in models:
            if ('seqlevel_featembeder_' in m_name) or \
               ('local_featemb_init_attn_' in m_name) or \
               ('local_featemb_mut_attn_' in m_name) or \
               ('global_featemb_init_attn_' in m_name) or \
               ('global_featemb_mut_attn_' in m_name) or \
               ('decoder_' in m_name):
                
                m_name_upd = "_".join(m_name.split('_')[:-1]) # remove the dataset suffix
                if base_model_suffix:
                    m_name_upd = f'{m_name_upd}_{base_model_suffix}'
            else:
                m_name_upd = m_name
            print('m_name:', m_name)
            print('m_name_upd:',m_name_upd)
            m.load_state_dict(torch.load(os.path.join(state_dict_dir, '{}.pkl'.format(m_name_upd)), map_location=device))

    # update models fdtype and move to device
    for m, m_name in models:
        m.type(fdtype).to(device)
    print('number of parameters of model:', get_num_trainable_params(models))
    # to freeze all layers except trainable_layernames
    freeze_layers(models, trainable_layernames)
    # get only parameters with require_grad=True
    # models_param = get_trainable_params(models)
    
    print('number of parameters of model after freezing:', get_num_trainable_params(models))

    c = 0
    for elm in get_trainable_params(models):
        c+=elm.numel()
    print(c)
    c = 0
    for elm in models_param:
        c+=elm.numel()
    print(c)
    print()

    if('train' in data_loaders):
        weight_decay = options.get('weight_decay', 1e-4)
        # see paper Cyclical Learning rates for Training Neural Networks for parameters' choice
        # `https://arxive.org/pdf/1506.01186.pdf`
        # pytorch version >1.1, scheduler should be called after optimizer
        # for cyclical lr scheduler, it should be called after each batch update
        num_iter = len(data_loaders['train'])  # num_train_samples/batch_size
        c_step_size = int(np.ceil(5*num_iter))  # this should be 2-10 times num_iter
        base_lr = 3e-4
        max_lr = 5*base_lr  # 3-5 times base_lr
        optimizer = torch.optim.Adam(models_param, weight_decay=weight_decay, lr=base_lr)
        cyc_scheduler = torch.optim.lr_scheduler.CyclicLR(optimizer, 
                                                          base_lr, 
                                                          max_lr, 
                                                          step_size_up=c_step_size,
                                                          mode='triangular', 
                                                          cycle_momentum=False)


    if ('validation' in data_loaders):
        m_state_dict_dir = create_directory(os.path.join(wrk_dir, 'model_statedict'))

    if(num_epochs > 1):
        fig_dir = create_directory(os.path.join(wrk_dir, 'figures'))

    # dump config dictionaries on disk
    config_dir = create_directory(os.path.join(wrk_dir, 'config'))
    ReaderWriter.dump_data(config, os.path.join(config_dir, 'mconfig.pkl'))
    ReaderWriter.dump_data(options, os.path.join(config_dir, 'exp_options.pkl'))
    # store attention weights for validation and test set
    seqid_fattnw_map = {dsettype: {} for dsettype in data_loaders if dsettype in {'test'}}
    mask_gen = MaskGenerator()

    for epoch in range(num_epochs):
        # print("-"*35)
        for dsettype in dsettypes:
            print("device: {} | experiment_desc: {} | run_num: {} | epoch: {} | dsettype: {} | pid: {}"
                  "".format(device, options.get('experiment_desc'), run_num, epoch, dsettype, pid))
            pred_score = []
            ref_score = []
            seqs_ids_lst = []
            dataset_ids_lst = []

            data_loader = data_loaders[dsettype]

            epoch_loss = 0.

            if(dsettype == 'train'):  # should be only for train
                for m, m_name in models:
                    m.train()
                requires_grad=True
            else:
                for m, m_name in models:
                    m.eval()
                requires_grad = False

            for i_batch, samples_batch_lst in enumerate(tqdm(data_loader)):
                loss_tensor = torch.zeros(len(samples_batch_lst)).type(fdtype).to(device)
                for i_data, samples_batch in enumerate(samples_batch_lst):
                    # print('batch num:', i_batch)

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
                    
                    
                    with torch.set_grad_enabled(dsettype == 'train'):
                        X_init_batch = init_annot_embed(X_init_nucl, X_init_proto, X_init_pbs, X_init_rt)
                        X_mut_batch = mut_annot_embed(X_mut_nucl, X_mut_pbs, X_mut_rt)
                        # print('X_init_batch.shape:', X_init_batch.shape)
                        # print('X_mut_batch.shape:',X_mut_batch.shape)
                        # print('x_init_len.shape', x_init_len.shape)
                        # print('x_mut_len.shape:', x_mut_len.shape)

                        # print(np.unique(x_init_len))
                        # print(np.unique(x_mut_len))
                        # (bsize,)
                        y_batch = y_val.type(fdtype).to(device)
                        
                        # masks (bsize, init_seqlen) or (bsize, mut_seqlen)
                        x_init_m = mask_gen.create_content_mask((X_init_batch.shape[0], X_init_batch.shape[1]), x_init_len)
                        x_mut_m =  mask_gen.create_content_mask((X_mut_batch.shape[0], X_mut_batch.shape[1]), x_mut_len)

                        __, z_init = init_encoder.forward_complete(X_init_batch, x_init_len, requires_grad=requires_grad)
                        __, z_mut =  mut_encoder.forward_complete(X_mut_batch, x_mut_len, requires_grad=requires_grad)
            
                        max_seg_len = z_init.shape[1]
                        init_mask = x_init_m[:,:max_seg_len].to(device)

                        # global attention
                        # s (bsize, embed_dim)
                        # assign each attention layer
                        if separate_attention_layers:
                            attn_module_indx = i_data
                        else:
                            attn_module_indx = 0

                        local_featemb_init_attn, __ = attn_modules_lst[attn_module_indx][0]
                        global_featemb_init_attn, __ = attn_modules_lst[attn_module_indx][1]
                        local_featemb_mut_attn, __ = attn_modules_lst[attn_module_indx][2]
                        global_featemb_mut_attn, __ = attn_modules_lst[attn_module_indx][3]
        

                        s_init_global, __ = global_featemb_init_attn(z_init, mask=init_mask)
                        # local attention
                        s_init_local, __ = local_featemb_init_attn(z_init, mask=X_init_rt[:,:max_seg_len])

                        max_seg_len = z_mut.shape[1]
                        mut_mask = x_mut_m[:,:max_seg_len].to(device)
                        s_mut_global, __ = global_featemb_mut_attn(z_mut, mask=mut_mask)
                        s_mut_local, __ = local_featemb_mut_attn(z_mut, mask=X_mut_rt[:,:max_seg_len])

                        seqfeat = seqlevel_featembeder(seqlevel_feat)
                        # y (bsize, 1)
                        # y_hat_logit, y_sigma = decoder(torch.cat([s_init_global, s_init_local, s_mut_global, s_mut_local, seqfeat], axis=-1))
                        
                        decoder, __ = decoder_lst[i_data]
                        y_hat_logit = decoder(torch.cat([s_init_global, s_init_local, s_mut_global, s_mut_local, seqfeat], axis=-1))

                        loss = loss_func(y_hat_logit, y_batch)
                        # print('loss.shape:',loss.shape)

                        loss_tensor[i_data] = loss.sum(axis=-1).mean()

                        pred_score.extend(torch.exp(y_hat_logit).tolist())
                        ref_score.extend(y_batch.tolist())
                        seqs_ids_lst.extend(list(b_seqs_id))
                        dataset_ids_lst.extend([datasets_name_lst[i_data]]*len(b_seqs_id))

                    # print('mean_loss:', loss)

                    # if loss_type != 'Logploss':
                    #     loss = loss_func(y_batch, y_hat_logit)
                    # else:
                    #     y_dist = torch.distributions.normal.Normal(y_hat_logit, y_sigma)
                    #     loss = loss_func(y_batch, y_dist)
                
                loss = loss_tensor.sum()

                if(dsettype == 'train'):
                    optimizer.zero_grad()
                    # print("computing loss")
                    # backward step (i.e. compute gradients)
                    # print('loss_tensor:', loss_tensor)
                    loss.backward()
                    # optimzer step -- update weights
                    optimizer.step()
                    # after each batch step the scheduler
                    cyc_scheduler.step()

                    # print(optimizer)

                    
                epoch_loss += loss.item()

            # end of epoch
            epoch_loss_avgbatch[dsettype].append(epoch_loss/len(data_loader))
            # get average modelscore across the datasets

            modelscore = perfmetric_report_multidata_cont(pred_score, 
                                                          ref_score, 
                                                          dataset_ids_lst,
                                                          seqs_ids_lst,
                                                          epoch_loss_avgbatch[dsettype][-1], 
                                                          epoch, flog_out[dsettype])
            
            perf = modelscore.spearman_corr
            # print('modelscore.spearman_corr:', modelscore.spearman_corr)
            # print('modelscore.pearson_corr:', modelscore.pearson_corr)
            # print('modelscore.spearman_lst:', modelscore.spearman_lst)
            # print('modelscore.pearson_lst:', modelscore.pearson_lst)
            # print('modelscore.modelscores_lst:', modelscore.modelscores_lst)

            if(perf > score_dict[dsettype].spearman_corr):
                # print(f'score_dict[{dsettype}].spearman_corr:', score_dict[dsettype].spearman_corr)

                score_dict[dsettype] = modelscore
                
                if(dsettype == 'validation'):
                    for m, m_name in models:
                        torch.save(m.state_dict(), os.path.join(m_state_dict_dir, '{}.pkl'.format(m_name)))
                    ReaderWriter.dump_data({'epoch':epoch+1}, os.path.join(m_state_dict_dir, 'best_epoch.pkl'))
                if dsettype in {'test', 'validation'}:
                    predictions_df = build_predictions_df(seqs_ids_lst, 
                                                          ref_score, 
                                                          pred_score, 
                                                          target_names,
                                                          dset_names=dataset_ids_lst)
                    predictions_path = os.path.join(wrk_dir, f'predictions_{dsettype}.csv')
                    predictions_df.to_csv(predictions_path)

    if(num_epochs > 1):
        for dsettype in epoch_loss_avgbatch:
            plot_xy(np.arange(num_epochs), 
                    epoch_loss_avgbatch[dsettype], 
                    'number of epochs', 
                    f'{loss_type} loss', 
                    'epoch batch average loss', 
                    dsettype,
                    fig_dir)
    # dump_scores
    dump_dict_content(score_dict, list(score_dict.keys()), 'score', wrk_dir)

    
def run_tune_cont_pe_RNN_kldiv(data_partition, dsettypes, config, options, wrk_dir, state_dict_dir=None, to_gpu=True, gpu_index=0):
    pid = "{}".format(os.getpid())  # process id description
    # get data loader config
    dataloader_config = config['dataloader_config']
    cld = construct_load_dataloaders(data_partition, dsettypes, dataloader_config, wrk_dir)
    # dictionaries by dsettypes
    data_loaders, epoch_loss_avgbatch, score_dict, flog_out = cld

    device = get_device(to_gpu, gpu_index)  # gpu device
    fdtype = options['fdtype']
    loss_type = options.get('loss_func', 'mse')

    if loss_type == 'MSEloss' or loss_type == 'RMSEloss':
        loss_func = nn.MSELoss(reduction='mean')
    elif loss_type == 'L1loss':
        loss_func = nn.L1Loss(reduction='mean')
    elif loss_type == 'Huberloss':
        loss_func = nn.SmoothL1Loss(reduction='mean')
    elif loss_type == 'KLDloss':
        loss_func = nn.KLDivLoss(reduction='none')


    num_epochs = options.get('num_epochs', 50)
    run_num = options.get('run_num')
    # parse config dict
    model_config = config['model_config']
    model_name = options['model_name']
    num_outcomes = options['num_outcomes']
    trainable_layernames = options['trainable_layernames']

    target_names = ['averageedited', 'averageunedited', 'averageindel']

    if(model_name == 'PE_RNN_kldiv'):
        annot_embed = options.get('annot_embed')
        assemb_opt = options.get('assemb_opt')
        seqlevel_featdim = options.get('seqlevel_featdim')

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

        local_featemb_init_attn = FeatureEmbAttention(z_dim)
        local_featemb_mut_attn = FeatureEmbAttention(z_dim)

        global_featemb_init_attn = FeatureEmbAttention(z_dim)
        global_featemb_mut_attn = FeatureEmbAttention(z_dim)

        seqlevel_featembeder = MLPEmbedder(inp_dim=seqlevel_featdim,
                                           embed_dim=z_dim,
                                           mlp_embed_factor=1,
                                           nonlin_func=model_config.nonlin_func,
                                           pdropout=model_config.p_dropout, 
                                           num_encoder_units=1)

        decoder  = MLPDecoderDistribution(5*z_dim,
                                        embed_dim=z_dim,
                                        outp_dim=num_outcomes,
                                        mlp_embed_factor=2,
                                        nonlin_func=model_config.nonlin_func, 
                                        pdropout=model_config.p_dropout, 
                                        num_encoder_units=1)

        # define optimizer and group parameters
        models = ((init_annot_embed, 'init_annot_embed'), 
                  (mut_annot_embed, 'mut_annot_embed'),
                  (init_encoder, 'init_encoder'),
                  (mut_encoder, 'mut_encoder'),
                  (local_featemb_init_attn, 'local_featemb_init_attn'),
                  (local_featemb_mut_attn, 'local_featemb_mut_attn'),
                  (global_featemb_init_attn, 'global_featemb_init_attn'),
                  (global_featemb_mut_attn, 'global_featemb_mut_attn'),
                  (seqlevel_featembeder, 'seqlevel_featembeder'),
                  (decoder, 'decoder'))

        models_param = []
        for m, m_name in models:
            if m_name in trainable_layernames:
                models_param.extend(list(m.parameters()))

            # init_params_(m)


    if(state_dict_dir):  # load state dictionary of saved models
        print('loading state dictionary from: ', state_dict_dir)
        for m, m_name in models:
            m.load_state_dict(torch.load(os.path.join(state_dict_dir, '{}.pkl'.format(m_name)), map_location=device))


    # update models fdtype and move to device
    print('casting to fdtype and moving to cuda!')
    for m, m_name in models:
        m.type(fdtype).to(device)

    # to freeze all layers except trainable_layernames
    freeze_layers(models, trainable_layernames)
    # get only parameters with require_grad=True
    # models_param = get_trainable_params(models)
    print(get_num_trainable_params(models))
    c = 0
    for elm in get_trainable_params(models):
        c+=elm.numel()
    print(c)
    c = 0
    for elm in models_param:
        c+=elm.numel()
    print(c)
    print()


    if('train' in data_loaders):
        weight_decay = options.get('weight_decay', 1e-4)
        # see paper Cyclical Learning rates for Training Neural Networks for parameters' choice
        # `https://arxive.org/pdf/1506.01186.pdf`
        # pytorch version >1.1, scheduler should be called after optimizer
        # for cyclical lr scheduler, it should be called after each batch update
        num_iter = len(data_loaders['train'])  # num_train_samples/batch_size
        c_step_size = int(np.ceil(5*num_iter))  # this should be 2-10 times num_iter
        base_lr = 3e-4
        max_lr = 5*base_lr  # 3-5 times base_lr
        optimizer = torch.optim.Adam(models_param, weight_decay=weight_decay, lr=base_lr)


        cyc_scheduler = torch.optim.lr_scheduler.CyclicLR(optimizer, 
                                                          base_lr, 
                                                          max_lr, 
                                                          step_size_up=c_step_size,
                                                          mode='triangular', 
                                                          cycle_momentum=False)


    if ('validation' in data_loaders):
        m_state_dict_dir = create_directory(os.path.join(wrk_dir, 'model_statedict'))

    if(num_epochs > 1):
        fig_dir = create_directory(os.path.join(wrk_dir, 'figures'))

    # dump config dictionaries on disk
    config_dir = create_directory(os.path.join(wrk_dir, 'config'))
    ReaderWriter.dump_data(config, os.path.join(config_dir, 'mconfig.pkl'))
    ReaderWriter.dump_data(options, os.path.join(config_dir, 'exp_options.pkl'))
    # store attention weights for validation and test set
    seqid_fattnw_map = {dsettype: {} for dsettype in data_loaders if dsettype in {'test'}}
    mask_gen = MaskGenerator()
    print('updated!!')
    for epoch in range(num_epochs):
        # print("-"*35)
        for dsettype in dsettypes:
            print("device: {} | experiment_desc: {} | run_num: {} | epoch: {} | dsettype: {} | pid: {}"
                  "".format(device, options.get('experiment_desc'), run_num, epoch, dsettype, pid))
            pred_score = []
            ref_score = []
            seqs_ids_lst = []
            data_loader = data_loaders[dsettype]

            epoch_loss = 0.

            if(dsettype == 'train'):  # should be only for train
                for m, m_name in models:
                    m.train()
                requires_grad=True
            else:
                for m, m_name in models:
                    m.eval()
                requires_grad = False

            for i_batch, samples_batch in enumerate(tqdm(data_loader)):
                # print('batch num:', i_batch)

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
                
                # zero gradient
                if(dsettype == 'train'):
                    optimizer.zero_grad()
                
                with torch.set_grad_enabled(dsettype == 'train'):
                    X_init_batch = init_annot_embed(X_init_nucl, X_init_proto, X_init_pbs, X_init_rt)
                    X_mut_batch = mut_annot_embed(X_mut_nucl, X_mut_pbs, X_mut_rt)
                    # print('X_init_batch.shape:', X_init_batch.shape)
                    # print('X_mut_batch.shape:',X_mut_batch.shape)
                    # print('x_init_len.shape', x_init_len.shape)
                    # print('x_mut_len.shape:', x_mut_len.shape)

                    # print(np.unique(x_init_len))
                    # print(np.unique(x_mut_len))
                    # (bsize,)
                    y_batch = y_val.type(fdtype).to(device)
                    
                    # masks (bsize, init_seqlen) or (bsize, mut_seqlen)
                    x_init_m = mask_gen.create_content_mask((X_init_batch.shape[0], X_init_batch.shape[1]), x_init_len)
                    x_mut_m =  mask_gen.create_content_mask((X_mut_batch.shape[0], X_mut_batch.shape[1]), x_mut_len)

                    __, z_init = init_encoder.forward_complete(X_init_batch, x_init_len, requires_grad=requires_grad)
                    __, z_mut =  mut_encoder.forward_complete(X_mut_batch, x_mut_len, requires_grad=requires_grad)
        
                    max_seg_len = z_init.shape[1]
                    init_mask = x_init_m[:,:max_seg_len].to(device)
                    # global attention
                    # s (bsize, embed_dim)
                    s_init_global, __ = global_featemb_init_attn(z_init, mask=init_mask)
                    # local attention
                    s_init_local, __ = local_featemb_init_attn(z_init, mask=X_init_rt[:,:max_seg_len])

                    max_seg_len = z_mut.shape[1]
                    mut_mask = x_mut_m[:,:max_seg_len].to(device)
                    s_mut_global, __ = global_featemb_mut_attn(z_mut, mask=mut_mask)
                    s_mut_local, __ = local_featemb_mut_attn(z_mut, mask=X_mut_rt[:,:max_seg_len])

                    seqfeat = seqlevel_featembeder(seqlevel_feat)
                    # y (bsize, 1)
                    # y_hat_logit, y_sigma = decoder(torch.cat([s_init_global, s_init_local, s_mut_global, s_mut_local, seqfeat], axis=-1))
                    
                    y_hat_logit = decoder(torch.cat([s_init_global, s_init_local, s_mut_global, s_mut_local, seqfeat], axis=-1))

                    loss = loss_func(y_hat_logit, y_batch)
                    # print(loss)
                    loss = loss.sum(axis=-1).mean()
                    # print('mean_loss:', loss)

                    # if loss_type != 'Logploss':
                    #     loss = loss_func(y_batch, y_hat_logit)
                    # else:
                    #     y_dist = torch.distributions.normal.Normal(y_hat_logit, y_sigma)
                    #     loss = loss_func(y_batch, y_dist)

                    if(dsettype == 'train'):
                        # print("computing loss")
                        # backward step (i.e. compute gradients)
                        loss.backward()
                        # optimzer step -- update weights
                        optimizer.step()
                        # after each batch step the scheduler
                        cyc_scheduler.step()

                        # print(optimizer)

                    pred_score.extend(torch.exp(y_hat_logit).tolist())
                    ref_score.extend(y_batch.tolist())
                    seqs_ids_lst.extend(list(b_seqs_id))

                    epoch_loss += loss.item()

            # end of epoch
            epoch_loss_avgbatch[dsettype].append(epoch_loss/len(data_loader))
            modelscore = perfmetric_report_cont(pred_score, ref_score, 
                                                epoch_loss_avgbatch[dsettype][-1], 
                                                epoch, flog_out[dsettype])
            perf = modelscore.spearman_corr
            if(perf > score_dict[dsettype].spearman_corr):
                score_dict[dsettype] = modelscore
                if(dsettype == 'validation'):
                    for m, m_name in models:
                        torch.save(m.state_dict(), os.path.join(m_state_dict_dir, '{}.pkl'.format(m_name)))
                        ReaderWriter.dump_data({'epoch':epoch+1}, os.path.join(m_state_dict_dir, 'best_epoch.pkl'))
                if dsettype in {'test', 'validation'}:
                    predictions_df = build_predictions_df(seqs_ids_lst, ref_score, pred_score, target_names)
                    predictions_path = os.path.join(wrk_dir, f'predictions_{dsettype}.csv')
                    predictions_df.to_csv(predictions_path)

    if(num_epochs > 1):
        for dsettype in epoch_loss_avgbatch:
            plot_xy(np.arange(num_epochs), 
                    epoch_loss_avgbatch[dsettype], 
                    'number of epochs', 
                    f'{loss_type} loss', 
                    'epoch batch average loss', 
                    dsettype,
                    fig_dir)
    # dump_scores
    dump_dict_content(score_dict, list(score_dict.keys()), 'score', wrk_dir)


def tune_trainval_run(datatensor_partitions, config_map, train_val_dir, state_dict_dir, run_gpu_map, num_epochs=20):
    dsettypes = ['train', 'validation']
    mconfig, options = config_map
    options['num_epochs'] = num_epochs  # override number of epochs using user specified value
    options['train_flag'] = True

    for run_num in datatensor_partitions:
        # update options run num to the current run
        options['run_num'] = run_num
        data_partition = datatensor_partitions[run_num]
        # tr_val_dir = create_directory(train_val_dir)
        path = os.path.join(train_val_dir, 'train_val', 'run_{}'.format(run_num))
        wrk_dir = create_directory(path)
        # print(wrk_dir)
        config_num = os.path.basename(train_val_dir)
        print(f'config_{config_num}')


        # state_dict_pth = os.path.join(train_dir, 'model_statedict')

        if options.get('model_name') == 'PE_RNN_kldiv':
            run_tune_cont_pe_RNN_kldiv(data_partition, 
                                        dsettypes, 
                                        mconfig, 
                                        options, 
                                        wrk_dir,
                                        state_dict_dir=state_dict_dir,
                                        to_gpu=True, 
                                        gpu_index=run_gpu_map[run_num])
        elif options.get('model_name') == 'PE_RNN_distribution_multidata':
            run_tune_cont_pe_RNN_distribution_multidata(data_partition,
                                                         dsettypes, 
                                                         mconfig, 
                                                         options, 
                                                         wrk_dir, 
                                                        state_dict_dir=state_dict_dir,
                                                        to_gpu=True, 
                                                        gpu_index=run_gpu_map[run_num])

def train_test_partition(datatensor_partition, config_map, tr_val_dir, run_gpu_map, queue):
    """To use this in multiprocessing module"""
    mconfig, options = config_map
    num_epochs = options['num_epochs']
    print('number of epochs:', num_epochs)
    # note: datatensor_partition and run_gpu_map are 1-entry dictionaries
    gpu_index = list(run_gpu_map.values())[0]
    partition_index = list(run_gpu_map.keys())[0]
    print('-- partition_index:', partition_index, 'gpu_index:', gpu_index, '--')
    train_val_run(datatensor_partition, config_map, tr_val_dir, run_gpu_map, num_epochs=num_epochs)
    test_run(datatensor_partition, config_map, tr_val_dir, tr_val_dir, run_gpu_map, num_epochs=1)
    queue.put(gpu_index)

def train_val_run(datatensor_partitions, config_map, train_val_dir, run_gpu_map, statedict_dir=None, num_epochs=20):
    dsettypes = ['train', 'validation']
    mconfig, options = config_map
    options['num_epochs'] = num_epochs  # override number of epochs using user specified value
    options['train_flag'] = True

    for run_num in datatensor_partitions:
        # update options run num to the current run
        options['run_num'] = run_num
        data_partition = datatensor_partitions[run_num]
        # tr_val_dir = create_directory(train_val_dir)
        path = os.path.join(train_val_dir, 'train_val', 'run_{}'.format(run_num))
        wrk_dir = create_directory(path)
        # print(wrk_dir)
        config_num = os.path.basename(train_val_dir)
        print(f'config_{config_num}')

        
        if statedict_dir is not None:
            state_dict_pth = os.path.join(statedict_dir, 
                                          'train_val', 
                                          'run_{}'.format(run_num),
                                          'model_statedict')
            
            print('state_dict_pth:', state_dict_pth)

        if options.get('model_name') == 'PE_RNN_distribution':
            run_cont_pe_RNN_distribution(data_partition, 
                               dsettypes, 
                               mconfig, 
                               options, 
                               wrk_dir,
                               state_dict_dir=None,
                               to_gpu=True, 
                               gpu_index=run_gpu_map[run_num])

def test_run(datatensor_partitions, config_map, train_val_dir, test_dir, run_gpu_map, num_epochs=1):
    dsettypes = ['test']
    mconfig, options = config_map
    options['num_epochs'] = num_epochs  # override number of epochs using user specified value
    options['train_flag'] = False
    for run_num in datatensor_partitions:
        # update options fold num to the current fold
        options['run_num'] = run_num
        data_partition = datatensor_partitions[run_num]
        train_dir = create_directory(os.path.join(train_val_dir, 'train_val', 'run_{}'.format(run_num)))
        if os.path.exists(train_dir):
            # load state_dict pth
            state_dict_pth = os.path.join(train_dir, 'model_statedict')
            path = os.path.join(test_dir, 'test', 'run_{}'.format(run_num))
            test_wrk_dir = create_directory(path)

                
            if options.get('model_name') == 'PE_RNN_distribution':
                run_cont_pe_RNN_distribution(data_partition, 
                                   dsettypes, 
                                   mconfig, 
                                   options, 
                                   test_wrk_dir,
                                   state_dict_dir=state_dict_pth, 
                                   to_gpu=True, 
                                   gpu_index=run_gpu_map[run_num])
                            
        else:
            print('WARNING: train dir not found: {}'.format(path))



def train_val_multidata_run(dtensor_partitions_multidata, config_map, train_val_dir, run_gpu_map, statedict_dir=None, num_epochs=20):
    dsettypes = ['train', 'validation']
    mconfig, options = config_map
    options['num_epochs'] = num_epochs  # override number of epochs using user specified value
    options['train_flag'] = True
    for run_num in dtensor_partitions_multidata:
        # update options run num to the current run
        options['run_num'] = run_num
        data_partition_lst = dtensor_partitions_multidata[run_num]
        # tr_val_dir = create_directory(train_val_dir)
        path = os.path.join(train_val_dir, 'train_val', 'run_{}'.format(run_num))
        wrk_dir = create_directory(path)
        # print(wrk_dir)
        config_num = os.path.basename(train_val_dir)
        print(f'config_{config_num}')


        if statedict_dir is not None:
            state_dict_pth = os.path.join(statedict_dir, 
                                          'train_val', 
                                          'run_{}'.format(run_num),
                                          'model_statedict')
            
            print('state_dict_pth:', state_dict_pth)
        else:
            state_dict_pth = None

        if options.get('model_name') == 'PE_RNN_distribution_multidata':
            run_cont_pe_RNN_distribution_multidata(data_partition_lst, 
                                                    dsettypes, 
                                                    mconfig, 
                                                    options, 
                                                    wrk_dir,
                                                    state_dict_dir=state_dict_pth,
                                                    to_gpu=True, 
                                                    gpu_index=run_gpu_map[run_num])



def test_multidata_run(dtensor_partitions_multidata, config_map, train_val_dir, test_dir, run_gpu_map, num_epochs=1):
    dsettypes = ['test']
    mconfig, options = config_map
    options['num_epochs'] = num_epochs  # override number of epochs using user specified value
    options['train_flag'] = False
    for run_num in dtensor_partitions_multidata:
        # update options fold num to the current fold
        options['run_num'] = run_num
        data_partition_lst = dtensor_partitions_multidata[run_num]
        train_dir = create_directory(os.path.join(train_val_dir, 'train_val', 'run_{}'.format(run_num)))
        if os.path.exists(train_dir):
            # load state_dict pth
            state_dict_pth = os.path.join(train_dir, 'model_statedict')
            path = os.path.join(test_dir, 'test', 'run_{}'.format(run_num))
            test_wrk_dir = create_directory(path)


            if options.get('model_name') == 'PE_RNN_distribution_multidata':
                run_cont_pe_RNN_distribution_multidata(data_partition_lst, 
                                                        dsettypes, 
                                                        mconfig, 
                                                        options, 
                                                        test_wrk_dir,
                                                        state_dict_dir=state_dict_pth, 
                                                        to_gpu=True, 
                                                        gpu_index=run_gpu_map[run_num])

        else:
            print('WARNING: train dir not found: {}'.format(path))