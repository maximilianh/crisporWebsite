'''
@author: orisenbazuru
'''
import torch
import torch.nn as nn
from torch.nn.utils.rnn import pack_padded_sequence, pad_packed_sequence

class RNN_Net(nn.Module):
    def __init__(self, 
                 input_dim, 
                 hidden_dim, 
                 z_dim, 
                 device,
                 num_hiddenlayers=1, 
                 bidirection= False, 
                 rnn_pdropout=0., 
                 rnn_class=nn.LSTM, 
                 nonlinear_func=nn.ReLU(),
                 fdtype = torch.float32):
        
        super().__init__()
        self.fdtype = fdtype
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.z_dim = z_dim
        self.num_hiddenlayers = num_hiddenlayers
        self.rnn_pdropout = rnn_pdropout
        self.device = device
        self.rnninput_dim = self.input_dim

        if num_hiddenlayers == 1:
            rnn_pdropout = 0
        self.rnn = rnn_class(self.rnninput_dim, 
                             hidden_dim, 
                             num_layers=num_hiddenlayers, 
                             dropout=rnn_pdropout, 
                             bidirectional=bidirection,
                             batch_first=True)
        if(bidirection):
            self.num_directions = 2
        else:
            self.num_directions = 1
   
        self.Wz = nn.Linear(self.num_directions*hidden_dim, self.z_dim)
        self.nonlinear_func = nonlinear_func    

        
    def init_hidden(self, batch_size, requires_grad=True):
        """initialize hidden vectors at t=0
        
        Args:
            batch_size: int, the size of the current evaluated batch
        """
        device = self.device
        # a hidden vector has the shape (num_layers*num_directions, batch, hidden_dim)
        h0=torch.zeros(self.num_hiddenlayers*self.num_directions, batch_size, self.hidden_dim).type(self.fdtype)
        h0.requires_grad=requires_grad
        h0 = h0.to(device)
        if(isinstance(self.rnn, nn.LSTM)):
            c0=torch.zeros(self.num_hiddenlayers*self.num_directions, batch_size, self.hidden_dim).type(self.fdtype)
            c0.requires_grad=requires_grad
            c0 = c0.to(device)
            hiddenvec = (h0,c0)
        else:
            hiddenvec = h0
        return(hiddenvec)
    
    def forward_tbptt(self, trunc_batch_seqs, hidden):
        # run truncated backprop
        trunc_rnn_out, hidden = self.rnn(trunc_batch_seqs, hidden)

        z_logit = self.nonlinear_func(self.Wz(trunc_rnn_out))
            
        return (hidden, z_logit)
    
    def detach_hiddenstate_(self, hidden):
        # check if hidden is not tuple # case of GRU or vanilla RNN
        if not isinstance(hidden, tuple):
            hidden.detach_()
        else: # case of LSTM
            for s in hidden:
                s.detach_()
    
    def forward_complete(self, batch_seqs, seqs_len, requires_grad=True):
        """ perform forward computation
        
            Args:
                batch_seqs: tensor, shape (batch, seqlen, input_dim)
                seqs_len: tensor, (batch,), comprising length of the sequences in the batch
        """

        # init hidden
        hidden = self.init_hidden(batch_seqs.size(0), requires_grad=requires_grad)
        # pack the batch
        packed_embeds = pack_padded_sequence(batch_seqs, seqs_len.cpu().numpy(), batch_first=True, enforce_sorted=False)
        packed_rnn_out, hidden = self.rnn(packed_embeds, hidden)

        # we need to unpack sequences
        unpacked_output, out_seqlen = pad_packed_sequence(packed_rnn_out, batch_first=True)
            
        z_logit = self.nonlinear_func(self.Wz(unpacked_output))
  
        return(hidden, z_logit)
    
    def forward(self, batch_seqs, seqs_len, requires_grad=True):
        return self.forward_complete(batch_seqs, seqs_len, requires_grad=requires_grad)