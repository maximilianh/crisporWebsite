SHELL = /bin/bash

All:
	make Fold;
	make partition;
	make Fold-smp;
	make partition-smp;

Fold: ../../exe/Fold Fold.out
Fold.out: 
	echo 'RR1664' > Fold.out; \
	time (Fold RR1664.seq RR1664.ct) 2>> Fold.out; \
	echo 'bac-ste' >> Fold.out; \
	time (Fold bac-ste.seq bac-ste.ct) 2>> Fold.out; \
	echo 'ivslus' >> Fold.out; \
	time (Fold ivslsu.seq ivslsu.ct) 2>> Fold.out; \
	echo 'ya5' >> Fold.out; \
	time (Fold ya5.seq ya5.ct) 2>> Fold.out; \
	echo 'ec16S' >> Fold.out; \
	time (Fold e.coli16s.25.4.seq ec16s.ct) 2>> Fold.out; \
	echo 'ec23S' >> Fold.out; \
	time (Fold ec23.seq bac-ste.ct) 2>> Fold.out; 

partition: ../../exe/partition partition.out
partition.out: 
	echo 'RR1664' > partition.out; \
	time ( ../../exe/./partition RR1664.seq RR1664.pfs ) 2>> partition.out; \
	echo 'bac-ste' >> partition.out; \
	time ( ../../exe/./partition bac-ste.seq bac-ste.pfs ) 2>> partition.out; \
	echo 'ivslsu' >> partition.out; \
	time ( ../../exe/./partition ivslsu.seq ivslsu.pfs ) 2>> partition.out; \
	echo 'ya5' >> partition.out; \
	time ( ../../exe/./partition ya5.seq ya5.pfs ) 2>> partition.out; \
	echo 'ec16S' >> partition.out; \
	time ( ../../exe/./partition e.coli16s.25.4.seq ec16s.pfs ) 2>> partition.out; \
	echo 'ec23S' >> partition.out; \
	time ( ../../exe/./partition ec23.seq ec23.pfs ) 2>> partition.out; 


Fold-smp: ../../exe/Fold-smp Fold-smp.out
Fold-smp.out: 
	echo 'RR1664' > Fold-smp.out; \
	time (Fold-smp RR1664.seq RR1664.ct) 2>> Fold-smp.out; \
	echo 'bac-ste' >> Fold-smp.out; \
	time (Fold-smp bac-ste.seq bac-ste.ct) 2>> Fold-smp.out; \
	echo 'ivslus' >> Fold-smp.out; \
	time (Fold-smp ivslsu.seq ivslsu.ct) 2>> Fold-smp.out; \
	echo 'ya5' >> Fold-smp.out; \
	time (Fold-smp ya5.seq ya5.ct) 2>> Fold-smp.out; \
	echo 'ec16S' >> Fold-smp.out; \
	time (Fold-smp e.coli16s.25.4.seq ec16s.ct) 2>> Fold-smp.out; \
	echo 'ec23S' >> Fold-smp.out; \
	time (Fold-smp ec23.seq bac-ste.ct) 2>> Fold-smp.out; 

partition-smp: ../../exe/partition-smp partition-smp.out
partition-smp.out: 
	echo 'RR1664' > partition-smp.out; \
	time (partition-smp RR1664.seq RR1664.pfs) 2>> partition-smp.out; \
	echo 'bac-ste' >> partition-smp.out; \
	time (partition-smp bac-ste.seq bac-ste.pfs) 2>> partition-smp.out; \
	echo 'ivslus' >> partition-smp.out; \
	time (partition-smp ivslsu.seq ivslsu.pfs) 2>> partition-smp.out; \
	echo 'ya5' >> partition-smp.out; \
	time (partition-smp ya5.seq ya5.pfs) 2>> partition-smp.out; \
	echo 'ec16S' >> partition-smp.out; \
	time (partition-smp e.coli16s.25.4.seq ec16s.pfs) 2>> partition-smp.out; \
	echo 'ec23S' >> partition-smp.out; \
	time (partition-smp ec23.seq bac-ste.pfs) 2>> partition-smp.out; 

clean:
	rm *.ct;rm *.pfs;

realclean:
	make clean;
	rm partition.out; rm Fold.out; rm Fold-smp.out; rm partition-smp.out
