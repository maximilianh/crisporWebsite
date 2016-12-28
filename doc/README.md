To test the program, first make sure that there is a directory "../genomes".
If it's not there, rename "genomes.sample" to "genomes":

    mv ../genomes.sample ../genomes

Then run this command:

    ../crispor.cgi sacCer3 sampleIn.fa sampleOut.mine.tsv -o sampleOutOfftargets.mine.tsv

The output file sampleOut.mine.tsv should be identical to sampleOut.tsv
sampleOutOfftargets.mine.tsv should be identical to sampleOutOfftargets.tsv

The files *.bowtie.tsv are the output files for the same command, just with the option --bowtie
in addition.

The file testInHg19.fa contains a sample for the hg19 genome, the output is in testOutHg19.tab 
and testOutHg19Offtargets.tab

    ../crispor.cgi hg19 testInHg19.fa testOutHg19.mine.tab -o testOutHg19Offtargets.mine.tab
