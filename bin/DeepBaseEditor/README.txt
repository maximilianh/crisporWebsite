1. System Requirements:
    Ubuntu   16.04
    Python   2.7.12
    Python Packages:
       numpy 1.14.5
       scipy 1.1.0

   Tensorflow and dependencies:
       Tensorflow  1.4.1
       CUDA       8.0.61
       cuDNN       6.0.21

2. Installation Guide (required time, <120 minutes):

- Operation System
    Ubuntu 16.04 download from https://www.ubuntu.com/download/desktop
    
 - Python and packages
    Download Python 2.7.12 tarball on https://www.python.org/downloads/release/python-2712/
    Unzip and install:
       tar -zxvf Python-2.7.12.tgz
       cd ./Python-2.7.12
       ./configure
       make

   Package Installation:
       pip install numpy==1.14.5
       pip install scipy==1.1.0

   Tensorflow Installation:
       (for GPU use)
       pip install tensorflow-gpu==1.4.1
       (for CPU only)
       pip install tensorflow==1.4.1


(for GPU use)

- CUDA Toolkit 8.0 
    wget -O cuda_8_linux.run https://developer.nvidia.com/compute/cuda/8.0/Prod2/local_installers/cuda_8.0.61_375.26_linux-run
    sudo chmod +x cuda_8_linux.run
    ./cuda_8_linux.run

- cuDNN 6.0.21
    Download CUDNN tarball on https://developer.nvidia.com/cudnn
    Unzip and install:
       tar -zxvf cudnn-8.0-linux-x64-v6.0.tgz

For more details, please refer to CUDA, CuDNN, and Tensorflow installation guide on Github:          
    https://gist.github.com/ksopyla/813a62d6afc4307755e5832a3b62f432


3. Demo Instructions (required time, <1 min):

(1) ABE_Efficiency

Input1: ./ABE_Efficiency_sample.txt        # List of Target Sequence(s)
    File format:
    Target number	30 bp target sequence (4 bp + 20 bp protospacer + PAM + 3 bp)			
	1	AACTGAAGGCTGAACAGCAGGGGTGGGAGA

Input2: ./ABE_Efficiency_Weight/ # Pre-trained Weight Files

Output: outputs/TEST_OUTPUT_for_test.xlsx
    Predicted activity score for sequence 1 :
    25.79517365

Run script:
    python ./TEST_ABE_Efficiency.py

Modification for personalized runs:

   <TEST_ABE_Efficiency.py>
    ## System Paths ##
    path                 = './'
    parameters           = {'0': 'ABE_Efficiency_sample.txt'}

   ## Run Parameters ##
    TEST_NUM_SET         = [0] # List can be expanded in case of multiple test parameters
    best_model_path_list = ['./ABE_Efficiency_Weight/']

ABE_Efficiency_Sample.txt can be replaced or modified to include target sequence of interest

 

(2) ABE_Proportion

Input1: ./ABE_Proportion_sample.txt        # List of Target Sequence(s)
    File format:
    Target number   30 bp target sequence (4 bp + 20 bp protospacer + PAM + 3 bp)
      1  		AACTGAAGGCTGAACAGCAGGGGTGGGAGA		
    
Input2: ./ABE_Proportion_Weight/ # Pre-trained Weight Files

Output: outputs/TEST_OUTPUT_for_test.xlsx
30 bp target sequence (4 bp + 20 bp protospacer + PAM + 3 bp)			Outcome seqeuence				Proportion
    AACTGAAGGCTGAACAGCAGGGGTGGGAGA				AACTGAGGGCTGAACAGCAGGGGTGGGAGA		0.068301663
    AACTGAAGGCTGAACAGCAGGGGTGGGAGA				AACTGAAGGCTGGACAGCAGGGGTGGGAGA		0.346454144
    AACTGAAGGCTGAACAGCAGGGGTGGGAGA				AACTGAGGGCTGGACAGCAGGGGTGGGAGA		0.006036451
    AACTGAAGGCTGAACAGCAGGGGTGGGAGA				AACTGAAGGCTGAGCAGCAGGGGTGGGAGA		0.433071852
    AACTGAAGGCTGAACAGCAGGGGTGGGAGA				AACTGAGGGCTGAGCAGCAGGGGTGGGAGA		0.023255052
    AACTGAAGGCTGAACAGCAGGGGTGGGAGA				AACTGAAGGCTGGGCAGCAGGGGTGGGAGA		0.11875882
    AACTGAAGGCTGAACAGCAGGGGTGGGAGA				AACTGAGGGCTGGGCAGCAGGGGTGGGAGA		0.000885288
	

Run script:
    python ./TEST_ABE_Proportion.py

Modification for personalized runs:

   <TEST_ABE_Proportion.py>
    ## System Paths ##
    path                 = './'
    parameters           = {'0': 'ABE_Proportion_sample.txt'}

   ## Run Parameters ##
    TEST_NUM_SET         = [0] # List can be expanded in case of multiple test parameters
    best_model_path_list = ['./ABE_Proportion_Weight/']

ABE_Proportion_sample.txt can be replaced or modified to include target sequence of interest

 

(3) CBE_Efficiency

Input1: ./CBE_Efficiency_sample.txt        # List of Target Sequence(s)
    File format:
    Target number   30 bp target sequence (4 bp + 20 bp protospacer + PAM + 3 bp)
  	    1   	TCAGGGCTGAACTAAAGCCTCCAGGGGGCC
    
Input2: ./CBE_Efficiency_Weight/ # Pre-trained Weight Files

Output: outputs/TEST_OUTPUT_for_test.xlsx
    Predicted activity score for sequence 1:
    4.856938124

Run script:
    python ./TEST_CBE_Efficiency.py

Modification for personalized runs:

   <TEST_CBE_Efficiency.py>
    ## System Paths ##
    path                 = './'
    parameters           = {'0': 'CBE_Efficiency_sample.txt'}

   ## Run Parameters ##
    TEST_NUM_SET         = [0] # List can be expanded in case of multiple test parameters
    best_model_path_list = ['./CBE_Efficiency_Weight/']

CBE_Efficiency_sample.txt can be replaced or modified to include target sequence of interest

 

(4) CBE_Efficiency_CA

Input1: ./CBE_Efficiency_CA_sample.txt        # List of Target Sequence(s)
    File format:
    Target number    30 bp target sequence (4 bp + 20 bp protospacer + PAM + 3 bp)		"Chromatin accessibility (1= DNase I hypersensitive sites, 0 = Dnase I non-sensitive sites)"
	1		TCAGGGCTGAACTAAAGCCTCCAGGGGGCC					1

Input2: ./CBE_Efficiency_CA_Weight/ # Pre-trained Weight Files

Output: outputs/TEST_OUTPUT_for_test.xlsx
    Predicted activity score for sequence 1:
    17.0684725

Run script:
    python ./TEST_CBE_Efficiency_CA.py

Modification for personalized runs:

   <TEST_CBE_Efficiency_CA.py>
    ## System Paths ##
    path                 = './'
    parameters           = {'0': 'CBE_Efficiency_CA_sample.txt'}

   ## Run Parameters ##
    TEST_NUM_SET         = [0] # List can be expanded in case of multiple test parameters
    best_model_path_list = ['./CBE_Efficiency_CA_Weight/']

CBE_Efficiency_CA_sample.txt can be replaced or modified to include target sequence of interest

 

(5) CBE_Proportion

Input1: ./CBE_Proportion_sample.txt        # List of Target Sequence(s)
    File format:
    Target number   30 bp target sequence (4 bp + 20 bp protospacer + PAM + 3 bp)
      1 		  TCAGGGCTGAACTAAAGCCTCCAGGGGGCC
  

Input2: ./CBE_Proportion_Weight/ # Pre-trained Weight Files

Output: outputs/TEST_OUTPUT_for_test.xlsx
    30 bp target sequence (4 bp + 20 bp protospacer + PAM + 3 bp)			Outcome seqeuence				Proportion   
	TCAGGGCTGAACTAAAGCCTCCAGGGGGCC			TCAGGGTTGAACTAAAGCCTCCAGGGGGCC		0.103429772
	TCAGGGCTGAACTAAAGCCTCCAGGGGGCC			TCAGGGCTGAATTAAAGCCTCCAGGGGGCC		0.846258879
	TCAGGGCTGAACTAAAGCCTCCAGGGGGCC			TCAGGGTTGAATTAAAGCCTCCAGGGGGCC		0.048627082


Run script:
    python ./TEST_CBE_Proportion.py

Modification for personalized runs:

   <TEST_CBE_Proportion.py>
    ## System Paths ##
    path                 = './'
    parameters           = {'0': 'CBE_Proportion_sample.txt'}

   ## Run Parameters ##
    TEST_NUM_SET         = [0] # List can be expanded in case of multiple test parameters
    best_model_path_list = ['./CBE_Proportion_Weight/']

CBE_Proportion_sample.txt can be replaced or modified to include target sequence of interest
