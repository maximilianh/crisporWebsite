To test the program, first make sure that there is a directory "../genomes".
If it's not there, rename "genomes.sample" to "genomes":

    mv ../genomes.sample ../genomes

Then run this command:

    ../crispor.cgi sacCer3 sampleIn.fa sampleOut.mine.tsv -o sampleOutOfftargets.mine.tsv

The output file sampleOut.mine.tsv should be identical to sampleOut.tsv
sampleOutOfftargets.mine.tsv should be identical to sampleOutOfftargets.tsv
