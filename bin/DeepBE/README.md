  


# DeepBE

Predicting base editing outcomes of 63 base editors

```diff
- If downloading models for use, please download them directly from the release note.
```

## System Requirements



```bash
    Python   3.6.13
    Python Packages:
    numpy 1.19.5
    pandas 1.3.0

    Tensorflow and dependencies:
    Tensorflow  2.6.2
    CUDA       11.2.0
    cuDNN       8.1.0
```
    
## Demo Instructions (required time, <1 min)

#### (1) DeepCas9variants
+ Predicting on-target efficiencies of 9 SpCas9 variants (SpCas9, VRQR variant, SpCas9-NG, SpCas9-NRRH, SpCas9-NRTH, SpCas9-NRCH, SpG, SpRY, Sc++)
+ PAM_variant feature: SpCas9 - 0, VRQR variant - 1, SpCas9-NG - 2, SpCas9-NRRH - 3, SpCas9-NRTH - 4, SpCas9-NRCH - 5, SpG - 6, SpRY - 7, Sc++ - 8

+ Input1: `./input_example.csv`  # List of Target Sequence(s)
    + File format:
```
    Target number	30 bp target sequence (4 bp + 20 bp protospacer + PAM + 3 bp), PAM_variant feature		
    GGATGACTACGCCTCTGCCTTagtAGGTCA, 2
```
+ Input2: `./PAM_variant_(*Cas9_variant_name)_model.h5` # Pre-trained Model Files

+ Output: ./prediction_result.xlsx
```
    target + PAM	feature	Prediction score
    GGATGACTACGCCTCTGCCTTagtAGGTCA	1	42.5612144470215
```
+ Run script:
```
    python ./PAM_variant_(*Cas9_variant_name).py
```

#### (2) DeepNG-BE
+ Predicting base editing outcomes of 7 base editors with SpCas9-NG
+ ABE: SpCas9-NG-ABE8e(V106W), SpCas9-NG-ABE8.17m+V106W
+ CBE: SpCas9-NG-YE1-BE4max, SpCas9-NG-SsAPOBEC3B
+ CGBE: SpCas9-NG-CGBE1, SpCas9-NG-miniCGBE1, SpCas9-NG-APOBEC-nCas9-Ung

+ Input1: `./input_example.csv`  # List of Target Sequence(s)
    + File format:
```
    Target number	30 bp target sequence (4 bp + 20 bp protospacer + PAM + 3 bp), PAM_variant feature		
    GGATGAACAACAAACTGCCTTAGTAGGTCA, 2
```
+ Input2: `./DeepNG-BE_(*base_editor_name)_model/` # Pre-trained Model Files

+ Output: ./prediction_result.xlsx
```
    target + PAM	edited output   PAM_variant   Prediction score
    GGATGAACAACAAACTGCCTTAGTAGGTCA	 GGATGAACAACAAgCTGCCTTAGTAGGTCA   2	0.00129620986990631
    GGATGAACAACAAACTGCCTTAGTAGGTCA	 GGATGAACAACAgACTGCCTTAGTAGGTCA   2	0.000881725805811584
    GGATGAACAACAAACTGCCTTAGTAGGTCA	 GGATGAACAACAggCTGCCTTAGTAGGTCA   2	0.0000012952218639839
    ... 
```
+ Run script:
```
    python ./DeepNG-BE_(*base_editor_name).py
```

#### (3) DeepBE
+ Predicting base editing outcomes of 63 base editors
+ SpCas9 variants: SpCas9, VRQR variant, SpCas9-NG, SpCas9-NRRH, SpCas9-NRTH, SpCas9-NRCH, SpG, SpRY, Sc++
+ ABE: ABE8e(V106W), ABE8.17m+V106W
+ CBE: YE1-BE4max, SsAPOBEC3B
+ CGBE: CGBE1, miniCGBE1, APOBEC-nCas9-Ung

+ Input1: `./input_example.csv`  # List of Target Sequence(s)
    + File format:
```
    Target number	30 bp target sequence (4 bp + 20 bp protospacer + PAM + 3 bp), PAM_variant feature		
    GGATGAACAACAAACTGCCTTAGTAGGTCA, 5
```
+ Input2: `./DeepBE_(*base_editor_name)_model/` # Pre-trained Model Files

+ Output: ./prediction_result.xlsx
```
Output: ./prediction_result.xlsx
    target + PAM	edited output   PAM_variant   Prediction score
    GGATGAACAACAAACTGCCTTAGTAGGTCA	 GGATGAACAACAAgCTGCCTTAGTAGGTCA   2	0.00635114079341292
    GGATGAACAACAAACTGCCTTAGTAGGTCA	 GGATGAACAACAgACTGCCTTAGTAGGTCA   2	0.0043202587403357
    GGATGAACAACAAACTGCCTTAGTAGGTCA	 GGATGAACAACAggCTGCCTTAGTAGGTCA   2	6.34630623608246E-06
    ... 
```
+ Run script:
```
    python ./DeepBE_(*base_editor_name).py
```
