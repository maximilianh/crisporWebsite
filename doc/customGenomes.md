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
```
