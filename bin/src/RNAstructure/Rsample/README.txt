Rsample is used to calculate the partition function for a sequence that has multiple conformations. 
To perform a complete Rsample calculation, the following steps should be performed: 

1. Run Rsample to produce a Partition Save File (PFS)
2. Run stochastic, with this PFS file as input, to produce a CT file with Boltzmann ensemble of 1,000 structures.
3. Read the CT file from step 2 (as a command line argument in Linux and MacOS) with the R script RsampleCluster.R
   (Located in the R-scripts folder).
   This program uses the algorithm by Ding and Lawrence to calculate optimal number of clusters 
   and their centroids. 
   In addition to R, it requires the installation of R package with clustering procedures 
   called fpc which can be done by typing install.packages("fpc") inside R.