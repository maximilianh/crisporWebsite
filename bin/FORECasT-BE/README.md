# FORECasT-BE

A tool for predicting guide efficacy for base editors. 

Input is a 20 nucletoide guide sequence (positions 1-20 where the PAM starts at position 21).

Output is a set of predictions for on-target base editing efficacy at positions 3-10 in the guide.

# Installation
First download the package:

```git clone https://github.com/ananth-pallaseni/FORECasT-BE.git```

Then run the following command inside the forecast-be directory:

```pip install -e .```

# Usage
```python
import forecast_be

target_seq = 'TCTGCTCAGCTCATGCCGAT' # 20nt spacer sequence 
editor_type = 'CBE' # or 'ABE'

# Predict the total fraction of edited reads for the target sequence
# If the `mean` and `std` arguments are None, then returns a z-score
# Input a mean and std to scale this into reael efficiency (good defaults are mean=0.5 & std=0.1)
guide_efficiency = forecast_be.predict_total(target_seq, editor=editor_type, mean=None, std=None)

# Predict the fraction of edited reads with the on-target substitituion at each position 
# Returns a list of predictions
positional_efficiency = forecast_be.predict_total(target_seq, editor=editor_type, mean=None, std=None)
```
