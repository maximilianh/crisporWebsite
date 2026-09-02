## global

- implementing SPiP (github.com/LBGC-CFB/SPiP) to predict disruption / creation of splice sites
- adding human cell lines reference genomes RPE1, K562, HAP1, HEK293T, cancer cell lines..
- choosing which DeepBE model to use (SpCas9 + 3 PAM variants with the best scores)
- adding base editing specific off-target scores (https://doi.org/10.1038/s41467-023-41004-3 - but the repo was deleted)
- adding the Jacquere library (https://doi.org/10.1016/j.xgen.2026.101190)
- staggered cut for eSpOT-ON (pam NGG-22) 
- in crisporAddGenome, replacement of gene models for NCBI / ENSEMBL genomes
- error 500 if an exon shorter than 23bp is entered in classic mode ?
- add custom PAM in all modes
- add custom tracks (in the test version so that the files can be reached)

## KO

- allow scanning for potential STOP codons with PAM variants :
    - add multipam mode in KO ?

## KI

- finish rescue mode
- add prime editing
- add double nicking strategy for KI with ssODN (Schubert et al. 2021) ~ In Progress
- adapt donor design rules for Cas12a (Schubert et al. 2021)
- si chargement de la page de résultats trop longue :
    - tableaux vides -> si clic affichage : requête -> AJAX -> remplissage du tableau
- simplify the display in showDonor () ?
- add CDS replacement ? https://doi.org/10.1038/s41467-023-42036-5
- Don't put dononr DNA sequence into the url (induces a size limit)

## notes 02/09/26

- Ajouter prime editing dans mode KO
- Biblio KI / Double nicking
- Renommer mode "Sequence not in referece genome" : ajouter alignement + message d'alerte + ajouter en mode classic
- RPE1 genome
- conserver tous les boutons en mode KI -> message si technique impossible
- vérifier sécurité IA
- captcha ?
- tuto
- manuel
