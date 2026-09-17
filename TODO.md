## global

- implementing SPiP (github.com/LBGC-CFB/SPiP) to predict disruption / creation of splice sites
    - for introduction of silent mutations on donor DNA
    - for taking into account silent bystander edits in precision editing mode
- adding human cell lines reference genomes RPE1, K562, HAP1, HEK293T, cancer cell lines..
- adding base editing specific off-target scores (https://doi.org/10.1038/s41467-023-41004-3 - but the repo was deleted)
- adding the Jacquere library (https://doi.org/10.1016/j.xgen.2026.101190)
- staggered cut for eSpOT-ON (pam NGG-22) 
- in crisporAddGenome, replacement of gene models for NCBI and ENSEMBL genomes
- add custom PAM in all modes
- Finish custom tracks
- add "custom base editor" mode (without scoring)

## KO 

- add removal of an out of frame exon

## KI

- si chargement de la page de résultats trop longue :
    - tableaux vides -> si clic affichage : requête -> AJAX -> remplissage du tableau
- add CDS replacement ? https://doi.org/10.1038/s41467-023-42036-5
- for hg19 / hg38, add a mode to fetch sequences from clinVar

## notes 02/09/26

- Ajouter prime editing dans mode KO
- Biblio KI / Double nicking
- Renommer mode "Sequence not in referece genome" : ajouter alignement + message d'alerte + ajouter en mode classic
- RPE1 genome
- vérifier sécurité -> injections + surcharge du serveur
- captcha ?
- tuto
- manuel
