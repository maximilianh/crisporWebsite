## global

- implementing SPiP (github.com/LBGC-CFB/SPiP) to predict disruption / creation of splice sites
- adding human cell lines reference genomes RPE1, K562, HAP1, HEK293T, cancer cell lines..
- choosing which DeepBE model to use (SpCas9 + 3 PAM variants with the best scores)
- adding base editing specific off-target scores
- adding the Jacquere library (https://doi.org/10.1016/j.xgen.2026.101190)
- staggered cut for eSpOT-ON (pam NGG-22) 
- in crisporAddGenome, replacement of gene models for NCBI / ENSEMBL genomes
- handle request when the query sequence is different from the genome
- error 500 if an exon shorted than 23bp is entered in classic mode ?

## KO

- extend KO with base editing to splice sites (DONE ?)
- allow scanning for potential STOP codons with PAM variants :
    - add multipam mode in KO ?

## KI

- finish rescue mode
- ajouter édition par prime editing
- renommer mode KI ("flexible editing ?")
- add double nicking strategy for KI with ssODN (Schubert et al. 2021)
