import torch
from torch import nn


class CELoss(nn.Module):
    def __init__(self, reduction='none'):
        super().__init__()
        self.reduction = reduction

    def forward(self, pred_log_prob, target_prob, weights=None):
        """
        Args:
            pred_log_prob: tensor, (bsize, num_outcomes), predicted log probability of outcomes
            target_prob: tensor, (bsize, num_outcomes), true probability of of outcomes
            weights: (optional), tensor, (bsize, num_outcomes), weights to be used for the different sequences
 
        """

        l = (-1*(target_prob*pred_log_prob))
        if weights:
            l = l*weights
        if self.reduction == 'sum':
            # (bsize,)
            return l.sum(axis=-1)
        elif self.reduction == 'mean':
            return l.sum(axis=-1).mean()
        else:
            # (bsize, num_outcomes)
            return l
        

### adapted from deepprime paper code
class BalancedMSELoss(nn.Module):

    def __init__(self, scale=True):
        super().__init__()

        #  'Correction_Deletion', 1, 3.81
        #  'Correction_Insertion', 2, 3.62
        #  'Correction_Replacement', 3, 2.17

        self.corrtype_mat_index = {'Correction_Deletion':0,
                                   'Correction_Insertion':1,
                                    'Correction_Replacement':2}
        self.corrtype_indx_wloss_map = {0: 0.6,
                                        1:0.7,
                                        2:1.}
        
        self.mse = nn.MSELoss()
        if scale:
            self.mse = ScaledMSELoss()
            print("Applying ScaledMSELoss")
        else:
            print("Applying MSELoss without scaling")

    def forward(self, pred, y,  correct_type_mat):
        pred = pred.view(-1, 1)
        y = y.view(-1, 1)
        # y = torch.log1p(actual[:, 0].view(-1, 1))
        total_loss = 0.
        for corr_type in ['Correction_Deletion', 'Correction_Insertion', 'Correction_Replacement']:
            ctype_indx = self.corrtype_mat_index[corr_type]
            ctype_factor = self.corrtype_indx_wloss_map[ctype_indx]
            ctype_cond = correct_type_mat[:, ctype_indx] == 1
            l = self.mse(pred[ctype_cond], y[ctype_cond]) * ctype_factor
            total_loss += l
        return total_loss


class ScaledMSELoss(nn.Module):

    def __init__(self):
        super().__init__()
        self.mseloss = nn.MSELoss(reduction='none')

    def forward(self, pred, y):
        """
        Args:
            pred: torch.tensor, shape (*, 1), prediction tensor from the model
            y: torch.tensor, shape (*, 1), based on torch.log1p(y) reference values
        """
        # Reciprocal of the square root of the original dataset distribution
        mu = torch.minimum(torch.exp(6 * (y-3)) + 1, torch.ones_like(y) * 5)

        return torch.mean(self.mseloss(pred, y) * mu)