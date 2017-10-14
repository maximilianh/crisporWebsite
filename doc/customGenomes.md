```
# requested by Yann Audic Rennes , email
# Xenopus laevis 9.1
wget ftp://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.1/1.8.3.2/XL_9.1_v1.8.3.2.primaryTranscripts.gff3.gz
wget ftp://ftp.xenbase.org/pub/Genomics/JGI/Xenla9.1/Xla.v91.repeatMasked.fa.gz
gunzip *.gz
sudo crisprAddGenome fasta Xla.v91.repeatMasked.fa --desc 'faXenLae91|Xenopus laevis|African Clawed Frog|Xenbase 9.1' --gff=XL_9.1_v1.8.3.2.primaryTranscripts.gff3
# Xenopus tropicalis 9.0
rm -f *.fa *.gff3
wget ftp://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.0/Xtropicalis.v9.repeatMasked.fa.gz
wget ftp://ftp.xenbase.org/pub/Genomics/JGI/Xentr9.0/Xtropicalisv9.0.Named.primaryTrs.gff3.gz
gunzip *.gz
sudo crisprAddGenome fasta *.fa --desc 'faXenTro90|Xenopus tropicalis|Western Clawed Frog|Xenbase 9.0' --gff=Xtropicalisv9.0.Named.primaryTrs.gff3 

# O. taurii Genbank v2, email from nigel.grimsley@obs-banyuls.fr
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/214/015/GCA_000214015.2_version_140606/GCA_000214015.2_version_140606_genomic.fna.gz
gunzip *.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/214/015/GCA_000214015.2_version_140606/GCA_000214015.2_version_140606_genomic.gff.gz
sudo crisprAddGenome fasta GCA_000214015.2_version_140606_genomic.fna --desc 'gbOstTau3|Ostreococcus tauri|O. tauri|Genbank 140606' --gff GCA_000214015.2_version_140606_genomic.gff 

# Bombus terrestris, email from igormmattos@gmail.com
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/255/GCF_000214255.1_Bter_1.0/GCF_000214255.1_Bter_1.0_genomic.fna.gz -q & 
wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/255/GCF_000214255.1_Bter_1.0/GCF_000214255.1_Bter_1.0_genomic.gff.gz &
sudo crisprAddGenome fasta *.fna --desc 'gbBomTer|Bombus terrestris|Buff-tailed bumblebee|Genbank Bter1.0' --gff *.gff 

# ferret, Scott Tyler, UIowa
sudo crisprAddGenome ucsc 9669

# Medicago, Pascal Ratet, CNRS/UPSUD Saclay, Thu Jan 26 19:02:07 CET 2017
wget http://www.medicagohapmap.org/downloads/R108/0.95/Medicago_R108_HM340_v0.95_assembly.fasta.gz
wget http://www.medicagohapmap.org/downloads/R108/0.95/Medicago_R108_HM340_v0.95_annotation.gff3.gz
sudo crisprAddGenome fasta *.fasta --desc 'faMedTru108|Medicago truncatula|barrelclover|Medicohapmap.org R108 0.95' 

# got genome by email from Emmanuelle
sed -i 's/|.*//g' *.fa  # remove | characters, not allowed in crispor
sudo crisprAddGenome fasta Slit-V2-gap-filled.fa --desc 'faSpoLit2|Spodoptera littoralis|African cotton leafworm|from Manue V2' --baseDir /data/www/crisporDev/crisporMax/genomes

# plasmoDB
wget http://plasmodb.org/common/downloads/Current_Release/PbergheiANKA/fasta/data/PlasmoDB-31_PbergheiANKA_Genome.fasta
wget http://plasmodb.org/common/downloads/Current_Release/PbergheiANKA/gff/data/PlasmoDB-31_PbergheiANKA.gff
sudo crisprAddGenome fasta PlasmoDB-31_PbergheiANKA_Genome.fasta --desc 'faPlaBer1|Plasmodium berghei ANKA|Plasmodium berghei ANKA|PlasmoDB 3.1' --baseDir /data/www/crispor/genomes --gff PlasmoDB-31_PbergheiANKA.gff 

wget http://plasmodb.org/common/downloads/Current_Release/PknowlesiH/fasta/data/PlasmoDB-31_PknowlesiH_Genome.fasta
wget http://plasmodb.org/common/downloads/Current_Release/PknowlesiH/gff/data/PlasmoDB-31_PknowlesiH.gff
sudo crisprAddGenome fasta PlasmoDB-31_PknowlesiH_Genome.fasta --desc 'faPlaKno1|Plasmodium knowlesi H|Plasmodium knowlesi H|PlasmoDB 3.1' --baseDir /data/www/crispor/genomes

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/188/465/GCA_001188465.1_ML2/GCA_001188465.1_ML2_genomic.gff.gz
gunzip *.gz
sudo crisprAddGenome fasta *.fna --desc 'gbMacLig2|Macrostomum lignano|Macrostomum lignano|Genbank GCA_001188465.1_ML2' --baseDir /data/www/crispor/genomes

# camelpox
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/839/105/GCF_000839105.1_ViralProj14156/GCF_000839105.1_ViralProj14156_genomic.gff.gz
sudo crisprAddGenome fasta *.fna --desc 'refCamelpox1|Camelpox Virus|Camelpox Virus|Refseq GCF_000839105.1' --baseDir /data/www/crispor/genomes

# manual download from dictybase.org
unzip dicty_gff3.zip 
cat usr/local/dicty/data/gff3/*.gff > dicty.gff
sudo crisprAddGenome fasta dicty_primary_genomic --desc 'faDicty1|Dictyostelium discoideum|Dictyostelium discoideum|Dictybase.org Apr 2017' --baseDir /data/www/crispor/genomes 

# Refseq goat genome - GFF files did not part correctly, so no gene models
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/704/415/GCF_001704415.1_ARS1/GCF_001704415.1_ARS1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/704/415/GCF_001704415.1_ARS1/GCF_001704415.1_ARS1_genomic.gff.gz
sudo crisprAddGenome fasta GCF_001704415.1_ARS1_genomic.fna --desc 'refCapHir1|Capra hircus|goat|RefSeq ARS1' --baseDir /data/www/crispor/genomes
wget ftp://ftp.kazusa.or.jp/pub/lotus/lotus_r3.0/Lj3.0_pseudomol.fna.gz
wget ftp://ftp.kazusa.or.jp/pub/lotus/lotus_r3.0/Lj3.0_gene_models.gff3.gz
sudo crisprAddGenome fasta Lj3.0_pseudomol.fna --desc 'faLotJap3|Lotus japonicus|Lotus japonicus|kazusa.or.jp V3' --baseDir /data/www/crispor/genomes

# got refseq files by email, but also at ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/558/735/GCA_001558735.1_ASM155873v1/GCA_001558735.1_ASM155873v1_assembly_report.txt
sudo crisprAddGenome fasta GCA_001558735.1_ASM155873v1_genomic.fna --desc 'refPanAgg1|Pantoea agglomerans|P. agg. FDAARGOS_160|Refseq GCF_001558735.1' --baseDir /data/www/crispor/genomes --gff GCA_001558735.1_ASM155873v1_genomic.gff 
```
sudo crisprAddGenome fasta *.fasta --desc 'tryBru42732|Trypanosoma brucei 427|T. brucei 427|TriTrypDb' --baseDir /data/www/crispor/genomes --gff *.gff 
#wget http://tritrypdb.org/common/downloads/Current_Release/TbruceiLister427/gff/data/TriTrypDB-32_TbruceiLister427.gff
#wget http://tritrypdb.org/common/downloads/Current_Release/TbruceiLister427/fasta/data/TriTrypDB-32_TbruceiLister427_Genome.fasta
#sudo crisprAddGenome fasta Pmamm_08Dec2013.fa --desc 'faPmammDec13|Phallusia mammillata|P. mammillata|C. Dantec, Dec 2013' --baseDir /data/www/crispor/genomes --gff Phmamm_transcripts_MTP2014_ANISEED.gff3 
sudo crisprAddGenome fasta Phfumi_2014_LengthScreen1000pb.txt --desc 'faPhFumi|Phallusia fumigata|P. fumigata|C. Dantec, 2014 LengthScreen' --baseDir /data/www/crispor/genomes --gff Phfumi_MTP2014_transcripts_filtered.gff3 
sudo crisprAddGenome fasta glossina-morsitans-yalescaffoldsgmory1.fa --desc 'faGloMorY16|Glossina morsitans|G. moritans|VectorBase Yale 1.6' --baseDir /data/www/crispor/genomes --gff glossina-morsitans-yalebasefeaturesgmory16.gff 
wget https://www.vectorbase.org/download/glossina-morsitans-yalebasefeaturesgmory16gff3gz
wget https://www.vectorbase.org/download/lutzomyia-longipalpis-jacobinabasefeaturesllonj13gff3gz
#wget https://www.vectorbase.org/download/lutzomyia-longipalpis-jacobinascaffoldsllonj1fagz
gunzip *.gz
#sudo crisprAddGenome fasta lutzomyia-longipalpis-jacobinascaffoldsllonj1fagz.fa --desc 'faLutLonLionJ1|Lutzomyia longipalpis|Sand fly|VectorBase LionJ1' --baseDir /data/www/crispor/genomes --gff lutzomyia-longipalpis-jacobinabasefeaturesllonj13gff3gz.gff 
mv anopheles-arabiensis-dongolabasefeaturesaarad16gff3gz anopheles-arabiensis-dongolabasefeaturesaarad16gff3gz.gff.gz
#wget https://www.vectorbase.org/download/anopheles-arabiensis-dongolabasefeaturesaarad16gff3gz
#wget https://www.vectorbase.org/download/anopheles-arabiensis-dongolascaffoldsaarad1fagz
mv anopheles-arabiensis-dongolascaffoldsaarad1fagz anopheles-arabiensis-dongolascaffoldsaarad1fagz.fa.gz
#sudo crisprAddGenome fasta CI_TypeB.fasta --desc 'cioIntB1|Ciona Intestinalis B|Sea squirt Type B|Agnes Roure' --baseDir /data/www/crispor/genomes 
sudo crisprAddGenome fasta sfru31.fa --desc 'faSpoFru1|Spodoptera frugiperda|Fall armyworm|E. Joly 3.1' --baseDir /data/www/crisporDev/crisporMax/genomes
wget https://bioinformatics.psb.ugent.be/gdb/ectocarpusV2/EctsiV2_genome.fasta.gz
cat *.gff3 > e.gff
sudo crisprAddGenome fasta *.fasta --desc 'faEctoSp2|Ectocarpus sp.|Ectocarpus sp.|V2 ugent.be' --baseDir /data/www/crispor/genomes --gff e.gff
wget https://www.vectorbase.org/download/biomphalaria-glabrata-bb02transcriptsbglab15fagz
mv biomphalaria-glabrata-bb02transcriptsbglab15fagz g.fa.gz
gunzip g.fa.gz 
wget https://www.vectorbase.org/download/biomphalaria-glabrata-bb02basefeaturesbglab15gff3gz -O annot.gff3.gz
wget https://www.vectorbase.org/download/biomphalaria-glabrata-bb02scaffoldsbglab1fagz -O genome.fa.gz
gunzip *.gz
sudo crisprAddGenome fasta *.fa --desc 'vecBioGla1|Biomphalaria glabrata|B. glabrata|VectorBase 1' --baseDir /data/www/crispor/genomes --gff annot.gff3
#sudo crisprAddGenome fasta *.fa --desc 'vecBioGla1|Biomphalaria glabrata|B. glabrata|VectorBase 1' --baseDir /data/www/crispor/genomes --gff annot.gff3
unzip Aspni7_download.zip 
mv Aspni7/* ./
gunzip *
sudo rm -rf /tmp2/pz12AspNig40
sudo crisprAddGenome fasta *.fa --desc 'pz12MarPol31|Marchantia polymorpha|Common liverwort|Phytozome 12, V3.1' --baseDir /data/www/crispor/genomes --gff Mpolymorpha_320_v3.1.gene.gff3  -h
sudo crisprAddGenome fasta Mpolymorpha_320_v3.0.fa --desc 'pz12MarPol31|Marchantia polymorpha|Common liverwort|Phytozome 12, V3.1' --baseDir /data/www/crispor/genomes --gff Mpolymorpha_320_v3.1.gene.gff3  -a 
