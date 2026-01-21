/*
 * A program that calculates the partition function for a strand of nucleic acids.
 * This strand of nucleic acids can be composed of either DNA or RNA.
 *
 * (c) 2009 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter
 */

#include "ProbScan_Interface.h"
#include <vector>
#include <string>
#include <algorithm>
///////////////////////////////////////////////////////////////////////////////
// Constructor.
//////////////////////////////////////////////////////////////////////////////
probscanInterface::probscanInterface() {

	

	// Initialize the calculation temperature.
	temperature = 310.15;

	// Initialize the nucleic acid type.
	alphabet = DT_RNA;

	//default threshold value; when searching across loops, they must exceed this threshold to be included in the output
	threshold = 0.01;

	//by default, used fixed width instead of scientific
	fixed = true;
}

///////////////////////////////////////////////////////////////////////////////
// Parse the command line arguments.
///////////////////////////////////////////////////////////////////////////////
bool probscanInterface::parse( int argc, char** argv ) {

	// Create the command line parser and build in its required parameters.
	ParseCommandLine* parser = new ParseCommandLine( "ProbScan" );
	parser->addParameterDescription( "input file", "The name of a file containing a partition function or input sequence." );

    //hairpin mode option
	vector<string> hairpinOptions;
	hairpinOptions.push_back( "-a" );
	hairpinOptions.push_back( "-A" );
	hairpinOptions.push_back( "--hairpin" );
	parser->addOptionFlagsNoParameters( hairpinOptions, "Print probabilities for all possible hairpin loops." );

    //internal loop mode option
	vector<string> iloopOptions;
	iloopOptions.push_back( "-i" );
	iloopOptions.push_back( "-I" );
	iloopOptions.push_back( "--internal" );
	parser->addOptionFlagsNoParameters( iloopOptions, "Print probabilities for all possible internal loops." );

    //bulge mode option
	vector<string> bulgeOptions;
	bulgeOptions.push_back( "-b" );
	bulgeOptions.push_back( "-B" );
	bulgeOptions.push_back( "--bulge" );
	parser->addOptionFlagsNoParameters( bulgeOptions, "Print probabilities for all possible bulge loops." );

    //stack mode option
	vector<string> helixOptions;
	helixOptions.push_back( "-e" );
	helixOptions.push_back( "-E" );
	helixOptions.push_back( "--helix" );
	parser->addOptionFlagsWithParameters( helixOptions, "Print probabilities for all possible helices with this number of base pair stacks. To get single base pair stacks, use -e 1." );

    //user specified loop mode option
	vector<string> pairOptions;
	pairOptions.push_back( "-p" );
	pairOptions.push_back( "-P" );
	pairOptions.push_back( "--pairs" );
	parser->addOptionFlagsWithParameters( pairOptions, std::string("Calculate probability for a user-specified loop. The loop must be provided as a set of pairs of nucleotide indices, where the nucs in the pair are delimited by dashes and each pair is delimited by a comma; eg \'-e 5-20\' will show the probability of a hairpin loop closed by a pair between nucleotides 5 and 20, and \'-e 10-120,15-70,75-110\' will give the probability of a three-way junction where the exiting helices are closed by pairs at 10-120, 15-70, and 75-110.") );

	vector<string> seqFileOptions;
	seqFileOptions.push_back( "-s" );
	seqFileOptions.push_back( "-S" );
	seqFileOptions.push_back( "--sequence" );
	parser->addOptionFlagsNoParameters( seqFileOptions, "Provide RNA from sequence file. Partition function will be calculated (may take a while); if you're going to query the same sequence repeatedly, you could save a lot of time by running from a partition function save file produced by the \'partition\' program." );

	vector<string> multibranchOptions;
	multibranchOptions.push_back( "-m" );
	multibranchOptions.push_back( "-M" );
	multibranchOptions.push_back( "--multibranch" );
	parser->addOptionFlagsWithParameters( multibranchOptions, "Provide a file with multibranch loops. These multibranch loops' probabilities will be checked." );

	vector<string> stemloopOptions;
	stemloopOptions.push_back("--stemloop");
	parser->addOptionFlagsWithParameters(stemloopOptions, std::string("Calculate probability for a hairpin stem-loop. The closing pair of the hairpin loop and the pair that closes the other end of the stem loop, where the nucs in the pair are delimited by dashes and each pair is delimited by a comma; eg \'--stemloop 100-105,98-107\' will show the probability of a hairpin loop closed by a pair between nucleotides 100 and 105, and with a helix with basepairs 98-107, 99-106, and 100-105."));


    // Add the DNA option.
    vector<string> dnaOptions;
    dnaOptions.push_back( "-d" );
    dnaOptions.push_back( "-D" );
    dnaOptions.push_back( "--DNA" );
    parser->addOptionFlagsNoParameters( dnaOptions, "Specify that the sequence is DNA, and DNA parameters are to be used. Default is to use RNA parameters." );

	// Add the scientific notation option.
	vector<string> fixedOptions;
	fixedOptions.push_back("--scientific");
	parser->addOptionFlagsNoParameters(fixedOptions, "Specify that output should be in scientific notation, rather than the default fixed notation.");

	// Add the Alphabet option.
	vector<string> alphabetOptions;
	alphabetOptions.push_back("--alphabet");
	parser->addOptionFlagsWithParameters(alphabetOptions, "Specify the name of a folding alphabet and associated nearest neighbor parameters. The alphabet is the prefix for the thermodynamic parameter files, e.g. \"rna\" for RNA parameters or \"dna\" for DNA parameters or a custom extended/modified alphabet. The thermodynamic parameters need to reside in the at the location indicated by environment variable DATAPATH. The default is \"rna\" (i.e. use RNA parameters). This option overrides the --DNA flag.");


	//threshold option
	vector<string> thresholdOptions;
	thresholdOptions.push_back("-t");
	thresholdOptions.push_back("-T");
	thresholdOptions.push_back("--threshold");
	parser->addOptionFlagsWithParameters(thresholdOptions, "Require the loops found in -a, -b, and -i modes to exceed this probability threshold to be included in the output.");

	// Parse the command line into pieces.
	parser->parseLine( argc, argv );

	// Get required parameters from the parser.
	if( !parser->isError() ) {
		inputFile = parser->getParameter( 1 );
	}
    //check the various modes
	if( !parser->isError() ) { hairpin = parser->contains( hairpinOptions ); }
	if( !parser->isError() ) {
        helix = parser->contains( helixOptions );
        parser->setOptionInteger(helixOptions, numstacks);
    }
	if( !parser->isError() ) { internal = parser->contains( iloopOptions ); }
	if( !parser->isError() ) { bulge = parser->contains( bulgeOptions ); }
	if( !parser->isError() ) {
        multibranch = parser->contains( multibranchOptions );
        loop_file = parser->getOptionString( multibranchOptions, true);
    }
	if( !parser->isError() ) {
        pairs = parser->contains( pairOptions );
        pairSpecification = parser->getOptionString( pairOptions, false);
    }
	if (!parser->isError()) {
		stemloop = parser->contains(stemloopOptions);
		pairstemloopSpecification = parser->getOptionString(stemloopOptions, false);
	}


	// Get the sequence file option.
	if( !parser->isError() ) { fromSequence = parser->contains( seqFileOptions ); }

	if (!parser->isError() && parser->contains(dnaOptions))
		alphabet = DT_DNA; // use DNA (unless overridden by alphabet)

	if (!parser->isError() && parser->contains(fixedOptions))
		fixed = false; // switch to scientific notation

	// Get the Alphabet option.
	if (!parser->isError() && parser->contains(alphabetOptions))
		alphabet = parser->getOptionString(alphabetOptions, false).c_str();

	// Get the Alphabet option.
	if (!parser->isError())
		parser->setOptionDouble(thresholdOptions, threshold);

	// Delete the parser and return whether the parser encountered an error.
	bool noError = ( parser->isError() == false );
	delete parser;
	return noError;
}


//functions to split an input string, for multibranch calculation IO
std::vector<std::string> &spl(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
//split string s in delimiter delim and return a vector of strings with the result
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    spl(s, delim, elems);
    return elems;
}
//convert string of the form "#-#" and returns a std::pair<int,int>
//eg "10-20" -> std::pair(10,20)
std::pair<int,int> stringtopair(const std::string s){
    std::vector<std::string> str_indices = split(s,'-');
//    cout<<s<<endl;
//    cout<<str_indices[0]<<endl;
//    cout<<str_indices[1]<<endl;
//    cout<<atoi(str_indices[0].c_str())<<endl;
    std::pair<int,int> p = std::make_pair(atoi(str_indices[0].c_str()),
                              atoi(str_indices[1].c_str()));
    return p;
}
//takes a string of form "#-#\t#-#\t ... #-#" and returns a vector of std::pair<int,int>
multibranch_loop_t stringtombl(const std::string s){
    std::vector<std::string> pairs = split(s,'\t');
    multibranch_loop_t mb;
    for(std::vector<std::string>::iterator it=pairs.begin()+1;it!=pairs.end();++it){
        mb.branches.push_back(stringtopair(*it));
    }
    return mb;
}

std::vector<std::pair<int,int> > processPairSpecification(std::string spec){
    std::vector<std::string> vs = split(spec,',');
    std::vector<std::pair<int,int> > vp;
    for(std::vector<std::string>::iterator it=vs.begin();it!=vs.end();++it){
        vp.push_back(stringtopair(*it));
    }
    std::sort(vp.begin(),vp.end());
    return vp;
}

bool forms_helix(std::vector<pair<int,int> > pairs)
{
    bool helix = true;
    for(size_t i=1;i<pairs.size();i++){
        std::pair<int,int> current = pairs[i];
        std::pair<int,int> last = pairs[i-1];
        if(!((current.first==last.first+1) && (current.second==last.second-1)))
            helix = false;
    }
    return helix;
}

bool inside(std::pair<int,int> outer,std::pair<int,int> inner)
{
    return ((outer.first<inner.first) && (outer.second>inner.second));
}

bool seperate(std::pair<int,int> a,std::pair<int,int> b){
    return (std::max(a.first,a.second) < std::min(b.first,b.second)) ||
        (std::max(b.first,b.second) < std::min(a.first,a.second));
}

bool forms_mbl(std::vector<pair<int,int> > pairs)
{
    bool mbl = true;
    std::pair<int,int> outer = pairs[0];
    for(size_t i=2;i<pairs.size();i++){
        std::pair<int,int> current = pairs[i];
        std::pair<int,int> last = pairs[i-1];
        if(!seperate(current,last)) mbl = false;
        if(!inside(outer,current)) mbl = false;
    }
    return mbl;
}

///////////////////////////////////////////////////////////////////////////////
// Run calculations.
///////////////////////////////////////////////////////////////////////////////
void probscanInterface::run() {
    
    if(!(hairpin||bulge||internal||multibranch||helix||pairs||stemloop)){
        std::cout<<"specify at least one type of loop, or ProbScan --help for usage message\n";
        return;
    }
    //main calculation
    try{
        ProbScan ps = ProbScan(inputFile.c_str(),fromSequence,alphabet.c_str());
        ErrorChecker<RNA> checker = ErrorChecker<RNA>( &ps );
        if(checker.isErrorStatus()) throw "error in initialization, did you try to initialize from a sequence without the -s flag?\n";
        if(hairpin){
            show_hairpins(ps.probability_of_all_hairpins(3,ps.GetSequenceLength()-2,threshold),fixed);
        }
        if(internal){
            show_internal_loops(ps.probability_of_all_internal_loops(threshold,std::string("internal")),fixed);
        }
        if(bulge){
            show_bulge_loops(ps.probability_of_all_internal_loops(threshold,std::string("bulge")),fixed);
        }
        if(helix){
            show_stacks(ps.probability_of_all_helices(threshold,numstacks),fixed);
        }
        if(multibranch){
            ifstream infile(loop_file.c_str());
            if (!infile.good()) throw "failed to open multibranch loop file\n";
            std::string line;
            while(getline(infile,line)){
                multibranch_loop_t mbl = stringtombl(line);
                mbl.probability = ps.probability_of_multibranch_loop(mbl);
                show_mbl(mbl,fixed);
            }
        }
        if(pairs){
            std::vector<std::pair<int,int> > p = processPairSpecification(pairSpecification);
            if(p.size()==0) throw "user specified loop is not well formed; use the form a-b,c-d\n";
            std::cout << "Probability of loop "<<pairSpecification<<": ";
            if(p.size()==1){
                std::cout<<ps.probability_of_hairpin(p[0].first,p[0].second);
            }
            else if(p.size()==2){
                int i = p[0].first;
                int j = p[0].second;
                int k = p[1].first;
                int l = p[1].second;
                if(i==k-1 && l==j-1){
                    std::cout<<ps.probability_of_helix(i,j,1);
                }
                else{
                    std::cout<<ps.probability_of_internal_loop(i,j,k,l);
                }
            }
            else{
                if (forms_helix(p)){
                    //process as helix
                    int i = p[0].first;
                    int j = p[0].second;
                    int k = p.back().first;
                    std::cout<<ps.probability_of_helix(i,j,k-i);
                }
                else if (forms_mbl(p)){
                    //make a multibranch loop
                    multibranch_loop_t mb = multibranch_loop(p[0].first,p[0].second);
                    for(std::vector<std::pair<int,int> >::iterator it = p.begin()+1;it!=p.end();++it){
                        mb.branches.push_back(*it);
                    }
                    std::cout<<ps.probability_of_multibranch_loop(mb);
                }
                else std::cout<<"...\nCould not interpret pair input, is it well formed?\nProvide pairs in the format i-j,i'-j', ... \nOne pair describes a hairpin,\ntwo pairs describes an internal loop or base pair stack,\nand three or more pairs describes a multibranch loop or helix.\nFor a helix, include ALL pairs that make up the helix\n";
            }
        }
		if (stemloop) {
			//User requested a stem-loop probability calculation
			std::vector<std::pair<int, int> > p = processPairSpecification(pairstemloopSpecification);
			if (p.size() != 2) throw "user specified stem-loop is not well formed; use the form a-b,c-d\n";
			else if (p[1].first- p[0].first!= p[0].second- p[1].second) throw "the number of pairs 5' and 3' to the hairpin loop are not equal\n";
			std::cout << "Probability of hairpin stem-loop " << pairstemloopSpecification << ": ";
			
			if (fixed) std::cout << std::fixed << std::setprecision(3) << ps.probability_of_stemloop(p[0].first, p[0].second, p[1].first, p[1].second);
			else std::cout << std::scientific << std::setprecision(3) << ps.probability_of_stemloop(p[0].first, p[0].second, p[1].first, p[1].second);
			

		}//end if stemloop
        cout<<std::endl;
        if(checker.isErrorStatus()) throw "error calculating loop probabilities\n";
    }
    catch(const char* oops){
        cout<<oops;
    }
}

///////////////////////////////////////////////////////////////////////////////
// Main method to run the program.
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] ) {

	probscanInterface* runner = new probscanInterface();
	bool parseable = runner->parse( argc, argv );
	if( parseable == true ) { runner->run(); }
	delete runner;
	return 0;
}
