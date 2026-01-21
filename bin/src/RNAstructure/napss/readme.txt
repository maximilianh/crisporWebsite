INTRODUCTION

The following information is for NAPSS version 1.0. (ca. 2008)
(Note that this information is outdated. The current NAPSS version includes several updates and modifications, 
including the ability to use triplet constraints. This file is kept for reference.)


The computer code contained herein is designed to enhance the accuracy of RNA secondary structure prediction by
incorporating experimental NMR data obtained for a particular RNA structure.  More specifically, these data 
describe helical regions of an RNA without knowledge of the exact nucleotides that comprise the basepairs within 
these helices.  Further information regarding this method can be found in the article "NMR-Assisted Prediction 
of RNA Secondary Structure: Identification of a Probable Pseudoknot in the Coding Region of an R2 
Retrotransposon" by Hart JM, Kennedy SD, Mathews DH, and Turner DH in the Journal of the American Chemical 
Society (2008) - in press at the time of this release.


LICENSING INFORMATION

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General 
Public License as published by the Free Software Foundation; either version 3 of the License, or (at your 
option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
GNU General Public License for more details.  You should have received a copy of the GNU General Public 
License along with this program (gpl.txt); if not, write to the Free Software Foundation, Inc., 59 Temple 
Place - Suite 330, Boston, MA 02111-1307, USA.


INSTALLATION

This software is a part of RNAstructure software package. Extract RNAstructure package to the target location
of your choice. The "napss" directory will contain the files necessary for executing the main NAPSS algorithm 
and some example data files. In order to successfully use the program an environment variable called DATAPATH 
must be set to the location of the data tables folder on the current machine. For example, in BASH, this is 
accomplished with: export DATAPATH=[directory in which RNAstructure resides]/RNAstructure/data_tables/ 
Execute a "make NAPSS" command in the RNAstructure directory to compile the NAPSS algorithm.
Execute a "make napss_dot" command in the RNAstructure directory  to compile the napss_dot algorithms.
napss_dot is the code to produce enhanced dotplots from an RNA sequence file, as well as some sample data files.


INPUT FILES

The NAPSS algorithm requires three text files as input: a sequence file from the program RNAstructure (also 
available from the Mathews Lab at http://rna.urmc.rochester.edu/software.html), a dotplot file created by the 
program "napss_dot" in the "RNAstructure/exe" directory of this release, and a user-created file containing the
helical walk constraints.  Additionally, it may be beneficial to create a fourth text file to specify the run-time 
configuration options for NAPSS.  Each of these files will be described in more detail below.  For simplicity,
it is recommended that users place each input file in the main NAPSS directory.

The sequence file is a text file that can be created in RNAstructure or a text editor.  It contains up to 
three comment lines at the beginning of the file that start with a semicolon.  This is followed by one line 
that contains the name of the sequence, and then as many lines as necessary to list the nucleotides 
(A/C/G/U) in the sequence, starting from the 5' end.  Immediately after the 3'-most nucleotide, the number 
"1" is appended to signal the end of the file.  Please see the file "bm.seq" in the napss directory.

The dotplot file is a tab-delineated text file created by the program in the "RNAstructure/exe" subdirectory of
this release.  To create this file, execute "napss_dot" program with the command "./napss_dot" followed by a space,
then the name of your sequence file, another space, and finally the desired name for the dotplot output.
For example, "./napss_dot bm.seq bm.dp". This procedure only needs to be performed once per unique RNA sequence;
multiple executions of the NAPSS algorithm on the same primary sequence can also use the same dotplot file. Please 
see the file "bm.dp" in the napss directory.

The constraint file should be generated with a plain text editor.  It describes the types of basepairs in 
each helical walk, and individual walks should be separated by carriage returns.  The current convention 
for numerically describing basepair types is as follows: "5" = AU, "6" = GC, and "7" = GU.  For example, 
the line "65666" would indicate a helical walk that consists of one GC pair followed by one AU pair 
followed by three more GC pairs.  Please see the file "bm.con" in the main directory for an example of 
this type of file.

The optional configuration file can also be generated with a plain text editor.  If this file is not 
specified at run-time, NAPSS will automatically enter an interactive mode to prompt the user for the 
information that this file would otherwise contain.  The configuration file consists of several lines of 
text with each line specifying one specific parameter in the format "PARAMETER=VALUE".  Please see the 
file "config.txt" in the main directory for an example of this type of file.
 
There are four mandatory parameters that a valid configuration file must contain: 
"inseq", the name of the sequence file 
"indotplot", the name of the dotplot file
"inconstraints", the name of the constraints file
"outct", the desired name for the connection-table file that describes the basepairs and calculated free 
energies of all the output secondary structures

Additionally, any combination of five optional parameters may also be specified: 
"maxtracebacks", the maximum number of refolded structures per constraint match combination [default is 
	100]
"percent", the maximum allowed percent difference from the lowest free energy structure in the dotplot 
	[default is 25]
"windowsize", a parameter describing how different suboptimal refoldings must be from each other 
	[default is 0] (a small window size allows very similar structures to be generated while a 
	larger window size requires them to be more different)
"cutoff", the maximum allowed percent difference from the lowest free energy structure in the final 
	output [default is 0, which directs NAPSS to output all structures (note that 0 is a flag
	that indicates "no cutoff," rather than a cutoff of 0)]
"outpairs", the desired name for the optional positions-paired output file for secondary structure 
	visualization with PseudoViewer [default is no output] (please see below for recommendations 
	on using this feature).


RUNNING NAPSS

Once all the input text files have been created, execution is accomplished by the command "./NAPSS",
optionally followed by a space and then by the name of the configuration text file.
For example, "./NAPSS config.txt"  If this file is not specified, NAPSS will automatically enter an
interactive mode to prompt the user for the necessary configuration parameters. Upon successful or 
unsuccessful completion, NAPSS will pause to allow Windows-based users to see the output before
terminating - press any character followed by Enter to return to the command prompt.


GENERAL RECOMMENDATIONS

The optional parameters for the configuration file have been assigned default values that correspond 
to the results reported in Hart JM, et al. JACS (2008), with the exception of "cutoff" which was 
assigned a value of 25.  In general the optimal values for each depend on the size of the RNA being 
studied as well as the length and uniqueness of the NMR-derived helical constraints.  Please refer 
to the specific discussion of these settings below.

The most critical parameter setting is "percent" - this, along with the size of the RNA, dictates 
how many potential basepairs from the dotplot are considered in the downstream calculations.  It is 
recommended that users start with a rather small value for this (5 or so) and slowly increase the 
value in subsequent calculations, particularly for structures that are larger than 100 nucleotides 
in length.  Larger "percent" values may be necessary to include poorly-predicted helices, especially 
in the case of certain pseudoknotted structures, but this can also exponentially increase the number 
of match combinations that NAPSS will have to refold.  Although this cannot at present be reduced to 
a mathematical relationship, a good rule of thumb is that a 100-nt RNA will take approximately one 
second per refolding on an average modern CPU.  NAPSS will update refolding and energy calculation 
progress every 100 structures so that users can arrive at a better performance estimate for their 
particular system.

"windowsize" should not be changed from the default value of zero without a compelling reason.  
This is because it is related only to the output of the refolding prediction, which tends to 
involve rather small substructures from unconstrained regions.  One potential exception to this 
would be if the overall structure is rather large and the number of constrained base pairs is 
comparatively small.

"maxtracebacks" governs the maximum number of structures that will be output from each refolding 
and can often be kept at a rather large number without any major detriment to running time.  It 
may be beneficial to decrease this setting if large regions of the secondary structure are 
unconstrained by experimental results, as this may decrease the number of relatively unstable 
structures while still outputting the most favorable candidates.

"cutoff" is an option for truncating the number of structures that will be output once the 
NAPSS algorithm has completed.  As such it does not have any effect on the running time of 
the algorithm but can potentially bring the number of results down to a more manageable size.  
One important caveat to note is that the free energy calculation that is currently incorporated 
in NAPSS cannot yield an accurate value for extremely complicated pseudoknotted structures (it 
is actually a variant of the Pseudoknot Energy Model in the NUPACK algorithm - for more 
information, please refer to Dirks RM, Pierce NA. J. Comput. Chem. 2003, 24, 1664-1677).  If 
such a structure is encountered by the energy-calculating subroutine, it breaks apart the 
offending substructure and adds a large penalty to that structure's predicted free energy.  If 
the consideration of such structures is desirable, "cutoff" should be set to zero so that all 
structures are output to the ct file.

The final optional parameter, "outpairs", can be used to specify a filename to which NAPSS will
 output all predicted secondary structures in a Positions-Paired format that can be read by the 
program PseudoViewer.  This program can be helpful for visualizng many pseudoknotted structures 
and is available for download from http://wilab.inha.ac.kr/pseudoviewer2/  The output file 
generated by this subroutine of NAPSS has many structures appended into one file with each 
structure separated from the next by a line of hyphens.  PseudoViewer can only read one 
structure at a time, so the desired structure should be copied and pasted into a new text 
file before attempting to visualize it in PseudoViewer.


BUG REPORTS

Please send a detailed email to James_Hart@urmc.rochester.edu with the subject line "NAPSS" 
to report any bugs that you find in this software.
