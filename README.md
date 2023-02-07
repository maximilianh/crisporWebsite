# CRISPOR - a CRISPR/Cas9 assistant 

CRISPOR predicts off-targets in the genome, ranks guides, highlights
problematic guides, designs primers and helps with cloning.  Try it on
http://crispor.org

CRISPOR uses BWA, a few tools from the UCSC Genome Browser (twoBitToFa, bedClip),
various R packages and a huge collection of external packages and source code files
from published articles, see the file crisporEffScores.py for the exact references or 
the tool tips when you mouse over the scores on the interactive website or the user's
manual http://crispor.org/manual/.

If you need to analyze hundreds of thousands of guides for a library, the tool FlashFry is
probably the better tool for you, see https://github.com/aaronmck/FlashFry. That being said, CRISPOR 
now has .bed input, so as long as you are not running on thousands of exons, it is probably
fast enough for most applications.

If you only need efficiency scores and no interactive website, try "python
crisporEffScores.py", it is a python module but also has a command line
interface that may be sufficient for programmers. There are also separate code files
for the Doench score and the MIT score that you can use for your projects, they 
don't have dependencies.

You can run crispor on the command line, or install it under an Apache webserver locally.
We also provide a virtual machine so you don't have to install anything yourself.

Note that usage is free only for academic or non-profit organisations, for commercial use
see license.txt. Licensing is managed by UCSC IP licensing, they have free demo licenses, 
various options adapted to your needs, by size of organisation. Licenses are
usually sold per-year and per-seat.

# Credits
* Jean-Paul Concordet had the original idea for CRISPOR. Without him the software would not exist.
* Alberto Stolfi for finding the N-SNP-bug
* Mark Diekhans for patching twoBitToFa and making it 100 times faster
* Many others! See the file changes.html for the full list of acknowledgements for every feature

# Software licenses

For the third-party software that you installed above, some of it is included for convenience under bin/:

* BWA is under GPL3
* libSVM: free and under copyright by Chih-Chung Chang and Chih-Jen Lin see http://www.csie.ntu.edu.tw/~cjlin/libsvm/COPYRIGHT
* svmlight: free for non-commercial use, see http://svmlight.joachims.org/
* SSC: no license specified
* primer3: GPL2.
* Fusi/Doench score: see LICENSE.txt, (c) by Microsoft Research

CRISPOR itself:

* The CRISPOR software, so the two Python files crispor.py and crisporEffScores.py, are released under a special license, see LICENSE.txt in this directory
