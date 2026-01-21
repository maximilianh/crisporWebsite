These are the RNAstructure text interfaces for Windows
for command line (DOS) use.

Note that the thermodynamic parameter files and HMM 
files, *.dat and *.dh files are required for operation.
You can either invoke a program from within these
directories as current working directory or set the 
DATAPATH environment variable to indicate the directory 
of their location.

An environment variable can be set on the command line 
with 
SET DATAPATH="C:\Documents and Settings\dhm\RNAstructure"

or (on Windows XP)
on the system properties dialog from the control panel 
under advanced and then clicking "Environment Variables."
Note that setting them here requires a reboot.



Most of these programs provide basic help if invoked
with no parameters and more extensive help when 
invoked with the -h flag.

Be sure to note the format of .seq files that are
input to RNAstructure.  There are several examples
included here.  These files require one or more comment 
lines that start with ";".  Then there is a sequence title 
line and finally the sequence, terminated with "1".

New in RNAstructure 5.01 is the ability to use FASTA 
sequence files.  Note the required format as shown 
in bmorivector.fasta.

Please report bugs to David Mathews:
David_Mathews@urmc.rochester.edu

