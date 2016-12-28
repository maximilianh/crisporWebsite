# CRISPOR - a CRISPR/Cas9 assistant 

CRISPOR predicts off-targets in the genome, ranks guides, highlights
problematic guides, designs primers and helps with cloning.  Try it on
http://crispr.org

CRISPOR uses BWA, a few tools from the UCSC Genome Browser (twoBitToFa, bedClip),
various R packages and a huge collection of external packages and source code files
from published articles, see the file crisporEffScores.py for the exact references.

If you only need efficiency scores and no interactive website, try "python
crisporEffScores.py", it is a python module but also has a command line
interface that may be sufficient for programmers.

# Installation of the command-line script:

Install BWA and a few required python modules:
    
    apt-get install bwa python-pip python-matplotlib
    sudo pip install biopython numpy scikit-learn==0.16.1 pandas twobitreader
    
or 
   
    yum install bwa python-pip
    sudo pip install biopython numpy scikit-learn==0.16.1 pandas matplotlib twobitreader
    
Install required R libraries:
   
    sudo Rscript -e 'install.packages(c("e1071"),  repos="http://cran.rstudio.com/")'
    sudo Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite(c("limma"));'

When you run crispor.py, it should show the usage message:
```
Usage: crispor.cgi [options] org fastaInFile guideOutFile 

Command line interface for the Crispor tool.

    org          = genome identifier, like hg19 or ensHumSap
    fastaInFile  = Fasta file
    guideOutFile = tab-sep file, one row per guide

    Use "noGenome" if you only want efficiency scoring (a LOT faster). This option 
    will use BWA only to match the sequence to the genome, extend it and obtain
    efficiency scores.

    If many guides have to be scored in batch: Add GGG to them to make them valid
    guides, separate these sequences by at least one "N" character and supply as a single
    fasta sequence, a few dozen to ~100 per file.
    

Options:
  -h, --help            show this help message and exit
  -d, --debug           show debug messages, do not delete temp directory
  -t, --test            run internal tests
  -p PAM, --pam=PAM     PAM-motif to use, default NGG. TTTN triggers special
                        Cpf1 behavior: no scores anymore + the PAM is assumed
                        to be 5' of the guide. Common PAMs are:
                        NGG,TTTN,NGA,NGCG,NNAGAA,NGGNG,NNGRRT,NNNNGMTT,NNNNACA
  -o OFFTARGETFNAME, --offtargets=OFFTARGETFNAME
                        write offtarget info to this filename
  -m MAXOCC, --maxOcc=MAXOCC
                        MAXOCC parameter, guides with more matches are
                        excluded
  --mm=MISMATCHES       maximum number of mismatches, default 4
  --bowtie              new: use bowtie as the aligner. Do not use. Bowtie
                        misses many off-targets.
  --skipAlign           do not align the input sequence. The on-target will be
                        a random match with 0 mismatches.
  --noEffScores         do not calculate the efficiency scores
  --minAltPamScore=MINALTPAMSCORE
                        minimum MIT off-target score for alternative PAMs, default
                        1.0
  --worker              Run as worker process: watches job queue and runs jobs
  --clear               clear the worker job table and exit
  -g GENOMEDIR, --genomeDir=GENOMEDIR
                        directory with genomes, default ./genomes
```
    
# Running the script as a CGI under Apache with the job queue

Make sure you can execute CGI scripts somewhere. Your Apache config (e.g. /etc/apache2/sites-enabled/000-default) should contain a section like this:

    <Directory "/var/www/html">
         AllowOverride All
         Options +ExecCGI (...)
         AddHandler cgi-script .cgi .pl .py

Also make sure you have the CGI module enabled:

    sudo a2enmod cgi
    sudo service apache2 restart

Clone the repo into such a directory:

    cd /var/www/html/
    git clone https://github.com/maximilianh/crisporWebsite
    
Use the sample E. coli genome for a start:

    mv genomes.sample genomes

Create a temp directory with the right permissions:
        
    mkdir temp
    chmod a+rw temp

By default, the jobs database is a SQlite file, /tmp/crisporJobs.db. The Apache
user has to be able to write to it so let us create it now:

    ./crispor.cgi --clear
    Worker queue now empty

Now start a single worker job. It will watch the job queue and process jobs:

    ./startWorkers.sh 1

Check that your worker is indeed running:
  
    cat log/worker1.log
    ps aux | grep crispor

Now try to access the script from a webbrowser, http://localhost/crispor.py and click "Submit"

# Adding a genome

Look into the "tools" directory [https://github.com/maximilianh/crisporWebsite/tree/master/tools], try the script crisprAddGenome.

The subdirectory usrLocalBin contains required tools for this script, you can copy them into /usr/local/bin of your machine, they are 64bit static linux binaries and should work on most current machines.

The script can auto-download genomes from Ensembl and UCSC or allows you to add your own custom genome in .fasta format. It does not handle gene models yet for custom genomes, email me if you need that, this step depends on the input file format of your genes.

E.g. to add the X. laevis genome:
    sudo crisprAddGenome fasta /tmp2/LAEVIS_7.1.repeatMasked.fa --desc 'xenBaseLaevis71|Xenopus laevis|X. laevis|Xenbase V7.1'

The four |-split values for the --desc option are: internalDatabaseName, scienticName, commonOrDisplayName, VersionOfAssembly

Make sure that internalDatabaseName does not include special characters, spaces etc.

# Thanks!
* Jean-Paul Concordet for numerous ideas on the user interface
* Alberto Stolfi for the finding the N-SNP-bug
* Mark Diekhans for patching twoBitToFa and making it 100 times faster

# Licenses

Included software:

* BWA is under GPL3
* libSVM: under copyright by Chih-Chung Chang and Chih-Jen Lin see http://www.csie.ntu.edu.tw/~cjlin/libsvm/COPYRIGHT
* svmlight: free for non-commercial use, see http://svmlight.joachims.org/
* SSC: no license specified
* primer3: GPL2.
* Fusi/Doench score: see LICENSE.txt, (c) by Microsoft Research
* crispor.py and crisporEffScores.py themselves are released under GPLv3, see LICENSE.txt
