## global

- implementing SPiP (github.com/LBGC-CFB/SPiP) to predict disruption / creation of splice sites
    - for introduction of silent mutations on donor DNA
    - for taking into account silent bystander edits in precision editing mode
- adding human cell lines reference genomes RPE1, K562, HAP1, HEK293T, cancer cell lines..
- adding the Jacquere library (https://doi.org/10.1016/j.xgen.2026.101190)
- adding base editing specific off-target scores (https://doi.org/10.1038/s41467-023-41004-3 - but the repo was deleted)
- staggered cut for eSpOT-ON (pam NGG-22) 
- in crisporAddGenome, replacement of gene models for NCBI and ENSEMBL genomes (currently ucsc only)
- add custom PAM in all modes
- Finish custom tracks
- add "custom base editor" mode (without scoring)
- create an agent to generate results from a prompt : https://www.nature.com/articles/s41586-026-11044-y 
- Add pegRNA design with OptiPrime : https://www.nature.com/articles/s41587-026-03261-7 - https://github.com/alvin-hsu/optiprime-src
- write the manual

## KO 

- add removal of an out of frame exon

## KI

- add CDS replacement ? https://doi.org/10.1038/s41467-023-42036-5
- for hg19 / hg38, add a mode to fetch sequences from clinVar
- optimize the loading of the results page ?
- simplify the display of the results page (above the sequence viewer)
- adapt donor synthesis constaints to the new twist standards
    - GC 50bp window 10-90%
    - 200bp direct repeats
    - 100bp hairpin forming repeats
    - 30bp homopolymers (all nucleotides)
    - 300-7000bp total length
- add an option to trim the donor DNA sequence to remove undesired features (e.g homopolymers) (only for double stranded donors)
- complete the list of tags and linkers
- if a gene annotation is available, codon-optimize the tag and linker sequences ?
- when re-designing pegRNAs witht silen bystanders, pass to regions to avoid (kozak + splice sites) in PRIDICT2 to save loading time
