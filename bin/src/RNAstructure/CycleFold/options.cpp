#include "options.h"
#include <iostream>
#include "../src/ParseCommandLine.h"
using std::string;
using std::vector;

options::options() : partition(false),
					 maxexpect(false),
					 turbo(false),
					 set(false),
					 iter(3),
					 gamma(0.6),
					 bigloops(false),
					 unpair(false)
{
}

options parse(int ac, char* av[]){
    options op = options();
    ParseCommandLine p("Options");

    p.addParameterDescription("sequence",
                              "The name of a file containing the sequence in fasta format.");

    vector<string> constraintOption = {"-c","-C","--constraintFile"};
    p.addOptionFlagsWithParameters(constraintOption, "Specify a constraint file to be applied (in constraint file format).");
    vector<string> ctOption = {"-ct","-CT","--constraintCt"};
    p.addOptionFlagsWithParameters(ctOption, "Specify a constraint file to be applied (in ct format).");
    vector<string> faConstraintOption = {"-fc","-FC","--fastaConstraints"};
    p.addOptionFlagsNoParameters(faConstraintOption, "Specify that the input fasta file contains secondary structure constraints (in dot-bracket format) to be applied to each structure. Default: off.");
    vector<string> unpairOp = {"-u","-U","--unpairingConstraints"};
    p.addOptionFlagsNoParameters(unpairOp, "Toggle whether restraints should be treated as unpairing constraints. Default: off.");
    vector<string> partitionOption = {"-p","-P","--partitionFunction"};
    p.addOptionFlagsNoParameters(partitionOption, "Specify that pair probabilities should be printed.");
    vector<string> maxExpect = {"-m","-M","--maxExpect"};
    p.addOptionFlagsNoParameters(maxExpect, "Specify that a MaxExpect calculation should be performed.");
    vector<string> turbo = {"-t","-T","--turbo"};
    p.addOptionFlagsNoParameters(turbo, "Specify that a TurboFold calculation should be performed.");
    vector<string> gamma = {"-g","-G","--gamma"};
    p.addOptionFlagsWithParameters(gamma, "Set gamma, the weighting parameter for extrinsic information inthe turbo calculation.");
    vector<string> iterations = {"-i","-I","--iterations"};
    p.addOptionFlagsWithParameters(iterations, "Set the number of iterations for the turbo calculation.");
    vector<string> big = {"-b","-B","--bigloops"};
    p.addOptionFlagsNoParameters(big, "Toggle whether large hairpins and internal loops are permitted in the structure. Default: off");
    vector<string> dotseqFormat = {"-s","-S","--seqFormat"};
    p.addOptionFlagsNoParameters(dotseqFormat, "Indicate that the input file is a SEQ formatted sequence (rather than a FASTA, which is the default).");

    p.parseLine(ac, av);
    op.input_files = {p.getParameter(1)};
    op.constraint_file = p.getOptionString(constraintOption, true);
    op.constraint_ct = p.getOptionString(ctOption, true);
    op.partition = p.contains(partitionOption);
    op.maxexpect = p.contains(maxExpect);
    op.turbo = p.contains(turbo);
    op.bigloops = p.contains(big);
    op.unpair = p.contains(unpairOp);
    op.fasta_constraints = p.contains(faConstraintOption);
	op.dotseq_format = p.contains(dotseqFormat);
    if(p.contains(gamma)) {
        op.gamma = atof(p.getOptionString(gamma, false).c_str());
      }
    else {
      op.gamma = 0.6;
    }
    if(p.contains(iterations)) {
        op.iter = atoi(p.getOptionString(iterations, false).c_str());
    }
    else {
      op.iter = 2;
    }
    if(!op.constraint_file.empty() && !op.constraint_ct.empty()){
      cout<<"warning: mutually exclusive options --constraintFile and --constraintCt chosen\n";
    }
	if (p.isHelp()) throw 0; // exit with return code 0 for the -h flag.
	if (p.isError()) throw "error parsing command line";
	const char* CYCLEFOLD_DATAPATH = getenv("CYCLEFOLD_DATAPATH");
	if(CYCLEFOLD_DATAPATH==NULL || string(CYCLEFOLD_DATAPATH)==string("")){
		throw "CYCLEFOLD_DATAPATH environment variable not set. Please set it to /path/to/RNAstructure/CycleFold/datafiles";
	}
	op.parameter_files = CYCLEFOLD_DATAPATH;
	if(op.parameter_files.back()!='/') op.parameter_files += string("/");
    return op;
}
