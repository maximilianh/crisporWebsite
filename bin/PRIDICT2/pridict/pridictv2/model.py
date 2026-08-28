import torch
from torch import nn

class MaskGenerator():
    def __init__(self):
        pass
    @classmethod
    def create_content_mask(clss, x_mask_shape, x_len):
        """
        Args:
            x_mask_shape: tuple, (bsize, max_seqlen)
            x_len: tensor, (bsize,), length of each sequence
        """
        x_mask = torch.ones(x_mask_shape)
        for bindx, tlen in enumerate(x_len):
            x_mask[bindx, tlen:] = 0
        return x_mask

class AnnotEmbeder_InitSeq(nn.Module):
    def __init__(self, embed_dim, annot_embed, assemb_opt='add'):
        super().__init__()
        self.num_nucl = 5 # nucleotide embeddings
        self.num_pbs_inidc = 2 # binary embedding
        self.num_annot_inidc = 3 # binary embedding

        self.assemb_opt = assemb_opt
        self.We = nn.Embedding(self.num_nucl+1, embed_dim, padding_idx=self.num_nucl)

        self.Wproto = nn.Embedding(self.num_annot_inidc+1, annot_embed, padding_idx=self.num_annot_inidc)
        self.Wpbs = nn.Embedding(self.num_pbs_inidc+1, annot_embed, padding_idx=self.num_pbs_inidc)
        self.Wrt = nn.Embedding(self.num_annot_inidc+1, annot_embed, padding_idx=self.num_annot_inidc)
    
    def forward(self, X_nucl, X_proto, X_pbs, X_rt):
        if self.assemb_opt == 'add':
            return self.We(X_nucl) + self.Wproto(X_proto) + self.Wpbs(X_pbs) + self.Wrt(X_rt)
        elif self.assemb_opt == 'stack':
            return torch.cat([self.We(X_nucl), self.Wproto(X_proto), self.Wpbs(X_pbs), self.Wrt(X_rt)], axis=-1)


class AnnotEmbeder_MutSeq(nn.Module):
    def __init__(self, embed_dim, annot_embed, assemb_opt='add'):
        super().__init__()

        self.num_nucl = 5 # nucleotide embeddings
        self.num_pbs_inidc = 2 # binary embedding
        self.num_annot_inidc = 3 # binary embedding

        self.assemb_opt = assemb_opt

        self.We = nn.Embedding(self.num_nucl+1, embed_dim, padding_idx=self.num_nucl)

        self.Wpbs = nn.Embedding(self.num_pbs_inidc+1, annot_embed, padding_idx=self.num_pbs_inidc)
        self.Wrt = nn.Embedding(self.num_annot_inidc+1, annot_embed, padding_idx=self.num_annot_inidc)
    
    def forward(self, X_nucl, X_pbs, X_rt):
        if self.assemb_opt == 'add':
            return self.We(X_nucl) + self.Wpbs(X_pbs) + self.Wrt(X_rt)
        elif self.assemb_opt == 'stack':
            return torch.cat([self.We(X_nucl), self.Wpbs(X_pbs), self.Wrt(X_rt)], axis=-1)


class SH_Attention(nn.Module):
    """ single head self-attention module
    """
    def __init__(self, input_size, embed_size):
        
        super().__init__()
        # define query, key and value transformation matrices
        # usually input_size is equal to embed_size
        self.embed_size = embed_size
        self.Wq = nn.Linear(input_size, self.embed_size, bias=False)
        self.Wk = nn.Linear(input_size, self.embed_size, bias=False)
        self.Wv = nn.Linear(input_size, self.embed_size, bias=False)
        self.softmax = nn.Softmax(dim=2) # normalized across feature dimension
        self.neginf = -1e6
    
    def forward(self, Xin_q, Xin_k, Xin_v, mask=None):
        """
        Args:
            Xin_q: query tensor, (batch, sequence length, input_size)
            Xin_k: key tensor, (batch, sequence length, input_size)
            Xin_v: value tensor, (batch, sequence length, input_size)
            mask: tensor, (batch, sequence length, sequence length) with 0/1 entries
                  (default None)
                  
        .. note:
            
            mask has to have at least one element in a row that is equal to one otherwise a uniform distribution
            will be genertaed when computing attn_w_normalized!
            
        """
        # print('---- SH layer ----')
        # print('Xin_q.shape', Xin_q.shape)
        # print('Xin_q.shape', Xin_k.shape)
        # print('Xin_v.shape', Xin_v.shape)

        # print('self.Wq:', self.Wq)
        # print('self.Wk:', self.Wk)
        # print('self.Wv:', self.Wv)

        X_q = self.Wq(Xin_q) # queries
        X_k = self.Wk(Xin_k) # keys
        X_v = self.Wv(Xin_v) # values
        
        # print('---- SH layer transform ----')
        # print('X_q.shape', X_q.shape)
        # print('X_k.shape', X_k.shape)
        # print('X_v.shape', X_v.shape)

        
        # scaled queries and keys by forth root 
        X_q_scaled = X_q / (self.embed_size ** (1/4))
        X_k_scaled = X_k / (self.embed_size ** (1/4))
        
        # (batch, sequence length, sequence length)
        attn_w = torch.bmm(X_q_scaled, X_k_scaled.transpose(1,2))
        # print('attn_w.shape:', attn_w.shape)
        # print()
         
        if mask is not None:
            # (batch, seqlen, seqlen)
            # if mask.dim() == 2: # assumption mask.shape = (seqlen, seqlen)
            #     mask = mask.unsqueeze(0) # add batch dimension
            # fill with neginf where mask == 0  
            attn_w = attn_w.masked_fill(mask == 0, self.neginf)
            # print('attn_w masked:\n', attn_w)

        attn_w_normalized = self.softmax(attn_w)
        # print('attn_w_normalized.shape', attn_w_normalized.shape)
        # print('attn_w_normalized masked:\n', attn_w_normalized)
        
        if mask is not None:
            # for cases where the mask is all 0 in a row
            attn_w_normalized = attn_w_normalized * mask
        
        # reweighted value vectors
        z = torch.bmm(attn_w_normalized, X_v)
        
        return z, attn_w_normalized


class FeatureEmbAttention(nn.Module):
    def __init__(self, input_dim):
        '''
        Args:
            input_dim: int, size of the input vector (i.e. feature vector)
        '''

        super().__init__()
        self.input_dim = input_dim
        # use this as query vector against the transformer outputs
        self.queryv = nn.Parameter(torch.randn(input_dim, dtype=torch.float32), requires_grad=True)
        self.softmax = nn.Softmax(dim=1) # normalized across seqlen
        self.neg_inf = -1e6

    def forward(self, X, mask=None):
        '''Performs forward computation
        Args:
            X: torch.Tensor, (bsize, seqlen, feature_dim), dtype=torch.float32
        '''

        X_scaled = X / (self.input_dim ** (1/4))
        queryv_scaled = self.queryv / (self.input_dim ** (1/4))
        # using  matmul to compute tensor vector multiplication
        
        # (bsize, seqlen)
        attn_w = X_scaled.matmul(queryv_scaled)


        if mask is not None:
            # (batch, seqlen)
            # fill with neginf where mask == 0  
            attn_w = attn_w.masked_fill(mask == 0, self.neg_inf)
            # print('attn_w masked:\n', attn_w)

        attn_w_normalized = self.softmax(attn_w)
        # print('attn_w_normalized.shape', attn_w_normalized.shape)
        # print('attn_w_normalized masked:\n', attn_w_normalized)
        
        if mask is not None:
            # for cases where the mask is all 0 in a row
            attn_w_normalized = attn_w_normalized * mask


        # reweighted value vectors (in this case reweighting the original input X)
        # unsqueeze attn_weights_norm to get (bsize, 1, seqlen)
        # perform batch multiplication with X that has shape (bsize, seqlen, feat_dim)
        # result will be (bsize, 1, feat_dim)
        # squeeze the result to obtain (bsize, feat_dim)
        z = attn_w_normalized.unsqueeze(1).bmm(X).squeeze(1)
        
        # returns (bsize, feat_dim), (bsize, seqlen)
        return z, attn_w_normalized


class MLPEmbedder(nn.Module):
    def __init__(self,
                 inp_dim,
                 embed_dim,
                 mlp_embed_factor=2,
                 nonlin_func=nn.ReLU(), 
                 pdropout=0.3, 
                 num_encoder_units=2):
        
        super().__init__()
        
        self.We = nn.Linear(inp_dim, embed_dim, bias=True)
        encunit_layers = [MLPBlock(embed_dim,
                                   embed_dim,
                                   mlp_embed_factor,
                                   nonlin_func, 
                                   pdropout)
                          for i in range(num_encoder_units)]

        self.encunit_pipeline = nn.Sequential(*encunit_layers)

    def forward(self, X):
        """
        Args:
            X: tensor, float32, (batch, embed_dim) representing x_target
        """

        X = self.We(X)
        out = self.encunit_pipeline(X)
        return out

class MLPBlock(nn.Module):
            
    def __init__(self,
                 input_dim,
                 embed_dim,
                 mlp_embed_factor,
                 nonlin_func, 
                 pdropout):
        
        super().__init__()
        
        assert input_dim == embed_dim

        self.layernorm_1 = nn.LayerNorm(embed_dim)

        self.MLP = nn.Sequential(
            nn.Linear(input_dim, embed_dim*mlp_embed_factor),
            nonlin_func,
            nn.Linear(embed_dim*mlp_embed_factor, embed_dim)
        )
        self.dropout = nn.Dropout(p=pdropout)

    def forward(self, X):
        """
        Args:
            X: input tensor, (batch, sequence length, input_dim)
        """
        o = self.MLP(X)
        o = self.layernorm_1(o + X)
        o = self.dropout(o)
        return o

class MLPDecoder(nn.Module):
    def __init__(self,
                 inp_dim,
                 embed_dim,
                 outp_dim,
                 mlp_embed_factor=2,
                 nonlin_func=nn.ReLU(), 
                 pdropout=0.3, 
                 num_encoder_units=2,
                 infer_sigma=False):
        
        super().__init__()

        self.infer_sigma = infer_sigma
        
        self.We = nn.Linear(inp_dim, embed_dim, bias=True)
        encunit_layers = [MLPBlock(embed_dim,
                                   embed_dim,
                                   mlp_embed_factor,
                                   nonlin_func, 
                                   pdropout)
                          for i in range(num_encoder_units)]

        self.encunit_pipeline = nn.Sequential(*encunit_layers)

        self.W_mu = nn.Linear(embed_dim, outp_dim)
        self.softplus = nn.Softplus()

        if self.infer_sigma:
            self.W_sigma = nn.Linear(embed_dim, outp_dim)

    def forward(self, X):
        """
        Args:
            X: tensor, float32, (batch, embed_dim) representing x_target
        """

        X = self.We(X)
        out = self.encunit_pipeline(X)
        
        mu = self.W_mu(out)
        if not self.infer_sigma:
            # return self.softplus(mu.view(-1))
            return self.softplus(mu)
        else:
            logsigma  = self.W_sigma(out)
            sigma = 0.1 + 0.9 * self.softplus(logsigma)
            return mu, sigma

class MLPDecoderDistribution(nn.Module):
    def __init__(self,
                 inp_dim,
                 embed_dim,
                 outp_dim,
                 mlp_embed_factor=2,
                 nonlin_func=nn.ReLU(), 
                 pdropout=0.3, 
                 num_encoder_units=2):
        
        super().__init__()
        
        self.We = nn.Linear(inp_dim, embed_dim, bias=True)
        encunit_layers = [MLPBlock(embed_dim,
                                   embed_dim,
                                   mlp_embed_factor,
                                   nonlin_func, 
                                   pdropout)
                          for i in range(num_encoder_units)]

        self.encunit_pipeline = nn.Sequential(*encunit_layers)

        self.W_mu = nn.Linear(embed_dim, outp_dim)
        self.log_softmax = nn.LogSoftmax(dim=-1)

    def forward(self, X):
        """
        Args:
            X: tensor, float32, (batch, embed_dim) representing x_target
        """

        X = self.We(X)
        out = self.encunit_pipeline(X)

        mu = self.W_mu(out)
        log_mu = self.log_softmax(mu)
        return log_mu

def init_params_(model):
    for p_name, p in model.named_parameters():
        param_dim = p.dim()
        if param_dim > 1: # weight matrices
            nn.init.xavier_normal_(p)
        elif param_dim == 1: # bias parameters
            if p_name.endswith('bias'):
                nn.init.uniform_(p, a=0.01, b=1.0)