# README #

### What is this repository for? ###
* This is a repository associated with the paper "Pairwise library screen systematically interrogates *Staphylococcus aureus* Cas9 specificity in human cells", by Tycko *et al.*, published in Nature Communications, 2018.
* We used the method to characterize the specificity of S. aureus Cas9, and here provide code to analyze the pairwise library screen sequencing data and model the specificity of SaCas9.

### Instructions for making predictions using existing SaCas9 specificity model


* To make predictions for a single guide and a single target sequence, use the following command. Both the guide sequence and the target sequenced should be entered as DNA and should not include the PAM.

```
./predict_activity_single.py guide_seq target_seq
```

### Instructions for training new model 

* If using new training data, the script `filter_dataset.py` can be used to filter out negative controls and guide-target pairs with low RBC counts. In the command below, `min_rbcs` is the minimum value of RBC counts that will be kept in downstream analysis. (To generate training data from a pairwise library screen, see below. The reads from the 2018 SaCas9 specificity pairwise library screen are available on SRA and the annotated library design is in the Supplementary Data of the paper.) 

```
./filter_dataset.py input_sheet.csv min_rbcs out_sheet.csv
```

* The next step is to encode the output from a previous step into a format usable by the training script.

```
./encode_full_matrix.py out_sheet.csv data_matrix.txt
```

* Finally, run the model on the file generated in the previous step. A single guide length must be specified in the `guide_length` argument. If the input sheet contains multiple guide lengths, only the specified one will be used. Derived parameter values will be deposited in `out_param_values.txt`.

```
Rscript train_rstan_model.R data_matrix.txt guide_length out_param_values.txt
```

### Instructions for analyzing a pairwise library screen

* Code for analyzing a new PWL is in the directory `pairwise-library-analysis`
* The script `PWL_library_features_annotate.py` is used to annotate the relationships between the pairwise guide and targets in each library member. It identifies the number, position, and type of mismatches and bulges. For an example of the output, see Supplementary Data.
* The script `PWL_screen_analysis.py` is used to compute the indel frequency and barcode representation for each library member, from an input of PEAR-assembled FASTQ files.
* The script `PWL_make_dataframe.py` is used to compute the Off: On-target ratios for each library member and output a pandas dataframe, from an input of the previously-computed indel frequencies. This dataframe was used to generate the figures in the paper.



### Who do I talk to? ###

* Repo admin:  
  Michael Dinsmore [michael.dinsmore@editasmed.com](mailto:michael.dinsmore@editasmed.com)  
    
	  
* Corresponding authors:  
  Patrick Hsu [patrick@salk.edu](mailto:patrick@salk.edu)  
  Christopher Wilson [christopher.wilson@editasmed.com](mailto:christopher.wilson@editasmed.com)