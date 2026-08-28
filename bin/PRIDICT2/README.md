# PRIDICT2.0: PRIme editing guide RNA preDICTion 2.0

![PRIDICT logo](dataset/PRIDICT2_logo.jpg)

## Index

1. [Overview](#1-overview)
2. [Complementary Models](#2-complementary-models)
3. [Additional Resources](#3-additional-resources)
4. [Contact](#4-contact)
5. [Citation](#5-citation)
6. [Getting Started](#6-getting-started)
   - 6.1 [Installation using Anaconda (Linux, Mac OS or WSL)](#61-installation-using-anaconda-linux-mac-os-or-wsl)
   - 6.2 [Running PRIDICT2.0 in 'single' mode](#62-running-pridict20-in-single-mode)
   - 6.3 [Running PRIDICT2.0 in 'batch' mode](#63-running-pridict20-in-batch-mode)
7. [Add-ons](#7-add-ons)
   - 7.1 [Prediction of pegRNAs with silent bystander edits](#71-prediction-of-pegrnas-with-silent-bystander-edits)
   - 7.2 [Prediction of pegRNAs with flexible insertion/deletion locations](#72-prediction-of-pegrnas-with-flexible-insertiondeletion-locations)

## 1. Overview

[PRIDICT2.0](https://rdcu.be/dLu0f) is an advanced version of the original [PRIDICT](https://rdcu.be/c3IM5) model designed for predicting the efficiency of prime editing guide RNAs. This repository allows you to run the model locally. For details on advancements over the original model, refer to our published study ([Mathis et al., Nature Biotechnology, 2024](https://rdcu.be/dLu0f)) and the initial [BioRxiv preprint](https://www.biorxiv.org/content/10.1101/2023.10.09.561414v1). For comprehensive step-by-step instructions, including practical tips for high-throughput screening, see our detailed protocol ([Mathis et al., Nature Protocols, 2025](https://rdcu.be/eAIEQ)).

## 2. Complementary Models

- **ePRIDICT**: This model focuses on the influence of local chromatin context (K562) on prime editing efficiencies and is designed to complement PRIDICT2.0. [Access GitHub Repository](https://github.com/Schwank-Lab/epridict)
- **DeepPrime**: This is a complementary model from [Yu et al. 2023](https://www.cell.com/cell/pdf/S0092-8674(23)00331-8.pdf) providing pegRNA design efficiency predictions for edit types <= 3 bp. [Access GitHub Repository](https://github.com/hkimlab/DeepPrime) 

## 3. Additional Resources

- **Protocol Paper**: For detailed step-by-step instructions on using PRIDICT2.0 and ePRIDICT, including practical tips & tricks for high-throughput screening, see our comprehensive protocol: [Mathis et al., Nature Protocols, 2025](https://rdcu.be/eAIEQ)
- **Supplementary Files**: [Access Here](https://github.com/Schwank-Lab/epridict/tree/supplementary_files)
- **Web Application**: For an online version of PRIDICT2.0, visit [our webapp](https://pridict.it/).

## 4. Contact

For questions or suggestions, please either:
- Email us at [nicolas.mathis@pharma.uzh.ch](mailto:nicolas.mathis@pharma.uzh.ch)
- Open a GitHub issue

## 5. Citation

If find our work useful for your research please cite:
- [Mathis et al., Nature Biotechnology, 2024](https://rdcu.be/dLu0f) (PRIDICT2.0):
```bibtex
@article{mathis2024machine,
  title={Machine learning prediction of prime editing efficiency across diverse chromatin contexts},
  author={Mathis, Nicolas and Allam, Ahmed and Tálas, András, and Kissling, Lucas and Benvenuto, Elena and Schmidheini, Lukas and Schep, Ruben and Damodharan, Tanav and Balázs, Zsolt and Janjuha, Sharan and others},
  journal={Nature Biotechnology},
  year={2024},
  publisher={Nature Publishing Group},
  doi={10.1038/s41587-024-02268-2}
}
```
- [Mathis & Allam et al., Nature Biotechnology, 2023](https://rdcu.be/c3IM5) (PRIDICT):
```bibtex
@article{mathis2023predicting,
  title={Predicting prime editing efficiency and product purity by deep learning},
  author={Mathis, Nicolas and Allam, Ahmed and Kissling, Lucas and Marquart, Kim Fabiano and Schmidheini, Lukas and Solari, Cristina and Balázs, Zsolt and Krauthammer, Michael and Schwank, Gerald},
  journal={Nature Biotechnology},
  volume={41},
  number={8},
  pages={1151--1159},
  year={2023},
  publisher={Nature Publishing Group},
  doi={10.1038/s41587-022-01613-7}
}
```
- [Mathis et al., Nature Protocols, 2025](https://rdcu.be/eAIEQ) (Protocol):
```bibtex
@article{mathis2025systematic,
  title={Systematic pegRNA design with PRIDICT2.0 and ePRIDICT for efficient prime editing},
  author={Mathis, Nicolas and Marquart, Kim Fabiano and Allam, Ahmed and Krauthammer, Michael and Schwank, Gerald},
  journal={Nature Protocols},
  year={2025},
  publisher={Nature Publishing Group},
  doi={10.1038/s41596-025-01244-7}
}
```





## 6. Getting Started

### 6.1 Installation using Anaconda (Linux, Mac OS or WSL) 🐍
Windows is only supported via [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install).

The easiest way to install and manage Python packages on various OS platforms is through [Anaconda](https://docs.anaconda.com/anaconda/install/). Once installed, any package (even if not available on Anaconda channel) could be installed using pip. 

* Install [Anaconda](https://docs.anaconda.com/anaconda/install/) or [miniconda](https://docs.anaconda.com/miniconda/). (conda 22.11 or newer)
* Start a terminal and run:
    ```shell
    # clone PRIDICT2.0 repository
    git clone https://github.com/uzh-dqbm-cmi/PRIDICT2.git
    # navigate into repository
    cd PRIDICT2
    # create conda environment and install dependencies for PRIDICT2 (only has to be done before first run/install)
    conda env create -f pridict2_repo.yml
        
    # activate the created environment
    conda activate pridict2

    # after activating environment, pytorch has to be installed separately here:
    pip install torch==2.0.1 --index-url https://download.pytorch.org/whl/cpu
  	
    # run desired PRIDICT2.0 command (single or batch mode, described below)
    python pridict2_pegRNA_design.py single --sequence-name seq1 --sequence "GCCTGGAGGTGTCTGGGTCCCTCCCCCACCCGACTACTTCACTCTCTGTCCTCTCTGCCCAGGAGCCCAGGATGTGCGAGTTCAAGTGGCTACGGCCGA(G/C)GTGCGAGGCCAGCTCGGGGGCACCGTGGAGCTGCCGTGCCACCTGCTGCCACCTGTTCCTGGACTGTACATCTCCCTGGTGACCTGGCAGCGCCCAGATGCACCTGCGAACCACCAGAATGTGGCCGC"
    # results are stored in 'predictions' folder
    ```

* `PRIDICT2.0` environment only has to be installed once. When already installed, follow the following commands to use `PRIDICT2.0` again:
    ```shell
    # open Terminal/Command Line
    # navigate into repository
    # activate the created environment
    conda activate pridict2
    # run desired PRIDICT2.0 command (single or batch mode, described below)
    python pridict2_pegRNA_design.py single --sequence-name seq1 --sequence "GCCTGGAGGTGTCTGGGTCCCTCCCCCACCCGACTACTTCACTCTCTGTCCTCTCTGCCCAGGAGCCCAGGATGTGCGAGTTCAAGTGGCTACGGCCGA(G/C)GTGCGAGGCCAGCTCGGGGGCACCGTGGAGCTGCCGTGCCACCTGCTGCCACCTGTTCCTGGACTGTACATCTCCCTGGTGACCTGGCAGCGCCCAGATGCACCTGCGAACCACCAGAATGTGGCCGC"
    # results are stored in 'predictions' folder
    ```

#

### 6.2 Running PRIDICT2.0 in 'single' mode:
  ####  Required:
  -  `--sequence-name`: name of the sequene (i.e. unique id for the sequence)
  -  `--sequence`: target sequence to edit in quotes (format: `"xxxxxxxxx(a/g)xxxxxxxxxx"`; minimum of 100 bases up and downstream of brackets are needed; put unchanged edit-flanking bases *outside* of brackets (e.g. xxxT(a/g)Cxxx instead of xxx(TAC/TGC)xxx)
  ####  Optional:
  -  `--output-dir`: output directory where results are dumped on disk (default: [`./predictions`](predictions); directory must already exist before running)
  -  `--use_5folds`: Use all 5-folds trained models. Default is to use fold-1 model
  -  `--nicking`: Additionally, design nicking guides for edit (PE3) with DeepSpCas9 prediction ([Kim et al. 2019](https://www.science.org/doi/10.1126/sciadv.aax9249)).
  -  `--ngsprimer`: Additionally, design high-throughput sequencing (NGS) primers for edit based on [Primer3](https://primer3.org/) design.

Example command:
```shell
python pridict2_pegRNA_design.py single --sequence-name seq1 --sequence "GCCTGGAGGTGTCTGGGTCCCTCCCCCACCCGACTACTTCACTCTCTGTCCTCTCTGCCCAGGAGCCCAGGATGTGCGAGTTCAAGTGGCTACGGCCGA(G/C)GTGCGAGGCCAGCTCGGGGGCACCGTGGAGCTGCCGTGCCACCTGCTGCCACCTGTTCCTGGACTGTACATCTCCCTGGTGACCTGGCAGCGCCCAGATGCACCTGCGAACCACCAGAATGTGGCCGC"
``` 
#

### 6.3 Running PRIDICT2.0 in 'batch' mode:
For instructions on running PRIDICT2.0 in batch mode on an HPC cluster with SLURM (e.g. with >1000 inputs), see our [HPC Batch Guide](./hpc_batch/README.md).
  ####  Required:
  -  `--input-fname`: input file name - name of csv file that has two columns [`editseq`, `sequence_name`]. See [`batch_template.csv`](input/batch_template.csv) in the [`./input`](input) folder
  ####  Optional:
  -  `--input-dir` : directory where the input csv file is found on disk
  -  `--output-dir`: directory on disk where to dump results (default: [`./predictions`](predictions))
  -  `--use_5folds`: Use all 5-folds trained models. Default is to use fold-1 model
  -  `--cores`: Number of batch input sequences that can be predicted simultaneously. Maximum 3 cores due to memory limitations. Default value 0 uses 3 cores if available.
  -  `--nicking`: Design nicking guides for edit (PE3) with DeepSpCas9 prediction ([Kim et al. 2019](https://www.science.org/doi/10.1126/sciadv.aax9249)).
  -  `--ngsprimer`: Additionally, design high-throughput sequencing (NGS) primers for edit based on [Primer3](https://primer3.org/) design.
  -  `--summarize`: Summarize the best scoring pegRNA(s) of each batch input in a summary file (saved in output folder in the format `date` +`time` + `HEK or K562` + `batch_summary.csv`). Choose either `HEK` or `K562` to get highest scoring pegRNAs based on desired PRIDICT2.0 score. (e.g., `--summarize K562`)
  -  `--summarize_number`: Define the number of top scoring pegRNAs to be included in summary file. Default is `3`.

Example command (including only required arguments):
```shell
 python pridict2_pegRNA_design.py batch --input-fname batch_template.csv
``` 
  #### Notes:
- A log file will be saved in [`log`](log) file folder (`date` +`time` + `batch filename (without file extension)` + `_batch_logfile.csv`). Check the column `log` to catch any inputs which contained errors. If everything worked well, it should read `Prediction successful!` in each row.
- To run PRIDICT2.0 from outside the repository folder, add path to script file (e.g. `/mnt/c/Users/pridictuser/github/PRIDICT2/pridict2_pegRNA_design.py`) and use `--input-dir` to define input folder (e.g. `/mnt/c/Users/pridictuser/github/PRIDICT2/input`)

--------------------------

## 7. Add-ons
### 7.1 Prediction of pegRNAs with silent bystander edits

`PRIDICT2.0`'s enhanced ability to predict multi-bp edits allows researchers to improve their outcomes by incorporating silent bystander edits into their targets. These edits can potentially boost editing efficiencies, particularly in MMR proficient contexts.

To facilitate this process, we provide a Jupyter Notebook [`notebook_silent_bystander_input.ipynb`](addons/silentbystander/notebook_silent_bystander_input.ipynb) in the [`addons/silentbystander`](addons/silentbystander/) folder. This notebook generates all possible `PRIDICT2.0` input sequences with silent bystanders up to 5bp up- and downstream of the edit.

#### Requirements:
- PRIDICT2.0 conda environment (includes necessary packages).
- PRIDICT input sequence with 150 bp context on both sides of the edit
- Input sequence must be in-frame (ORF_start = 0) if the edit and its context are within an exon
- If not in an exon, still assume it is "in frame" and also keep the `in frame value` as `yes` if you do batch input file
  --> to get all possible mutations for an intron/non-genic region, enter your input sequence in all ORF variants (3 for the fw strand; 3 for the rv strand) and run the batch mode.

#### Usage overview:
1. Use the notebook or corresponding command line script to create an input batch file for `PRIDICT2.0` prediction with silent bystanders.
2. Supports 1bp replacements and multi-bp replacements (not insertions/deletions).
3. Input format: `PRIDICT2.0` format with 150bp flanking bases on both sides.
4. Choose between single (single mutation/edit) or batch (multiple inputs) functions.
5. *Run `PRIDICT2.0` in batch mode in a terminal outside of notebook with the generated input sequences to obtain efficiency predictions.
(Optional: run PRIDICT2.0 from within notebook, see commented out example)*
6. Summarize individual bystander predictions to one prediction file (best prediction of each bystander variant)

#
### 7.2 Prediction of pegRNAs with flexible insertion/deletion locations

If the edit location is flexible (e.g., inserting a stop codon at the best site), multiple `PRIDICT2.0` predictions need to be performed to identify the most efficient option. To facilitate this process, we provide a Jupyter Notebook [`notebook_flexible_mutations.ipynb`](addons/flexible_mutations/notebook_flexible_mutations.ipynb) in the [`addons/flexible_mutations`](addons/flexible_mutations) folder. This notebook automates the generation of `PRIDICT2.0` input sequences for flexible mutations.

#### Requirements:
- PRIDICT2.0 conda environment (includes necessary packages).
- Target sequence with the possible insertion/deletion region enclosed in square brackets (`[ ]`; **not `()`**), including min. **100 bp** of context on both sides.
Example:
```shell
TGCCTGGAGGTGTCTGGGTCCCTCCCCCACCCGACTACTTCACTCTCTGTCCTCTCTGCCCAGGAGCCCAGGATGTGCGAGTTCAAGTGGCTACGGCCGA[CTGTCCTCTCTGCCCAGG]GTGCGAGGCCAGCTCGGGGGCACCGTGGAGCTGCCGTGCCACCTGCTGCCACCTGTTCCTGGACTGTACATCTCCCTGGTGACCTGGCAGCGCCCAGATG
```



For **insertions**:
- Set `"insertion"` as the `edit_type`.
- Define the inserted bases using IUPAC base codes (ATGC, N, R, etc.).
- Define insertion frequency (default = `1`; choose `3` for in-frame insertions).
- For in-frame insertions: Adjust the target sequence so that the **ORF starts at the brackets** (e.g., add 2bp at the beginning to reach 102bp, which is divisible by 3).

For **deletions**:
- Set `"deletion"` as the `edit_type`.
- Define the deletion length.
- Define deletion frequency (default = `1`; choose `3` for in-frame deletions).
- For in-frame deletions: Adjust the target sequence so that the **ORF starts at the brackets** (e.g., add 2bp at the beginning to reach 102bp, which is divisible by 3).

#### Usage overview:
1. Use this notebook to create an input batch file for `PRIDICT2.0` predictions with flexible mutations (insertions or deletions).
2. Input format: Target sequence with square brackets (`[ ]`) defining the region where insertion or deletion options should be created.
3. Choose between single (one condition) or batch (multiple conditions via `.csv` input; example file under [`addons/flexible_mutations/input/input_flexible_mutations_testfile.csv`](addons/flexible_mutations/input/input_flexible_mutations_testfile.csv) processing.
4. *Run `PRIDICT2.0` in batch mode in a terminal outside of notebook with the generated input sequences to obtain efficiency predictions.
(Optional: run PRIDICT2.0 from within notebook, see commented out example)*
5. Summarize predictions of all flexible mutation options into a single file by selecting the best predicted pegRNA for each variant.

---