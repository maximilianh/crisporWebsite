## global

- implementing SPiP (github.com/LBGC-CFB/SPiP) to predict disruption / creation of splice sites
- adding human cell lines reference genomes RPE1, K562, HAP1, HEK293T, cancer cell lines..
- choosing which DeepBE model to use (SpCas9 + 3 PAM variants with the best scores)
- adding base editing specific off-target scores
- adding the Jacquere library (https://doi.org/10.1016/j.xgen.2026.101190)
- staggered cut for eSpOT-ON (pam NGG-22) 
- in crisporAddGenome, replacement of gene models for NCBI / ENSEMBL genomes
- handle request when the query sequence is different from the genome
- error 500 if an exon shorter than 23bp is entered in classic mode ?
- add custom PAM mode in all modes

## KO

- extend KO with base editing to splice sites (DONE ?)
- allow scanning for potential STOP codons with PAM variants :
    - add multipam mode in KO ?

## KI

- finish rescue mode
- add prime editing
- add double nicking strategy for KI with ssODN (Schubert et al. 2021) ~ In Progress
- adapt donor design rules for Cas12a (Schubert et al. 2021)
- si chargement de la page de résultats trop longue :
    - tableaux vides -> si clic affichage : requête -> AJAX -> remplissage du tableau
