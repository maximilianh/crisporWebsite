/*
 * This is the command-line interface to Dynalign, 
 * written by David Mathews; Copyright 2002, 2003, 2004, 2005, 2006
 *
 * Contributors:
 *  Chris Connett and Andrew Yohn, 2006
 *  Josh Keegan, 2006
 *  Arif Harmanci, 2006
 *
 * Dynalign is described in:
 * Mathews & Turner, JMB, 317:191-203 (2002).
 *
 *
 * 
 *----------------------------------------------------------------
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdlib>

#include "../src/configfile.h"
#include "../src/defines.h"
#include "../src/dynalign.h"
#include "../src/dynalignheap.h"
#include "../src/outputconstraints.h"
#include "../src/platform.h"
#include "../src/rna_library.h"
#include "../src/structure.h"
#include "../src/TProgressDialog.h"
#include "../src/phmm/phmm_aln.h"
#include "../src/phmm/structure/structure_object.h"
#include "../RNA_class/RNA.h"
#include "../src/ErrorChecker.h"

//TOLERANCE is the maximum deviation from 310.15 K before which the enthalpy parameters 
//are read from disk to adjust the free energyy changes from 310.15 K.
#define TOLERANCE 0.01 

//#ifdef COMPILE_SMP  RMW 2016-06-11 - the interface for observingtextprogressbar has been changed.
//#include "../src/observingtextprogressbar.h"
//#endif

using namespace std;

// Zane: to get seq filename to pass to sanity check
#ifdef CHECK_ARRAY
string seq1;
string seq2;
#endif

bool openseq(structure& ct, const char*const seq) {
	int result=ct.openseqx(seq);
    if (result!=0)
		cerr<<"ERROR: Could not open sequence file "<<seq<<": "<<RNA::GetErrorMessage(result)<<ct.GetErrorDetails() << endl;
	return result==0;
}

//The main entry point for Dynalign using a text interface:
int main(int argc, char* argv[]) {
  string inseq1;
  string inseq2;
  string outct;
  string outct2;
  string aout;
  string alphabet;
  int imaxseparation;
#ifdef DYNALIGN_II
  double slope;
  double intercept;
  int max_elongation;
#else
#endif
  double fgap;
  int maxpairs=-1;

#ifdef DYNALIGN_II  
  int islope;
  int iintercept;
#else
#endif
  int igapincrease;
  int maxtrace;
  int percent;
  int singlefold_subopt_percent,maximumpairingdistance1,maximumpairingdistance2;
  int bpwin;
  int awin;
#ifdef DYNALIGN_II
#else
  bool insert;
#endif
  string savefile;
  string constraint1;
  string constraint2;
  string constrainta;
  string shape1;
  string shape2;
  bool dsv_templated, ct_templated;
  string dsvtemplatename;
	float maxdsvchange;
  int numProcessors = 1;
  bool **allowed_alignments,optimalonly;
  bool local;
  FILE *check;//a c FILE for checking that specified files exist
  
	int i,query,index;
	short **align;

	double shapeslope1,shapeintercept1,shapeslope2,shapeintercept2;

	double Temperature;
	bool DNA;


	TProgressDialog progress;
  
	bool constrained;
	short **forcealign;

	//#define cdbg cout << "DEBUG " << __LINE__ << ": "
  //cdbg << "TOP" << endl;
  
  // Check if we have a command-line parameter; if not, fall back to
  // interactive mode.
	if( argc == 2 ) {

    ConfigFile config(argv[1]);
    
    // Check for mandatory config file settings.
    bool valid_config =
      config.contains("inseq1") &&
      config.contains("inseq2") &&
      config.contains("outct") &&
      config.contains("outct2") &&
      config.contains("aout");
      //(config.contains("imaxseparation")||config.contains("align_threshold")) &&
      //config.contains("fgap") &&
      //config.contains("maxtrace") &&
      //config.contains("percent") &&
      //config.contains("singlefold_subopt_percent") &&
      //config.contains("bpwin") &&
      //config.contains("awin") &&
      //config.contains("insert") &&
      //config.contains("dsv_templated");

#ifdef COMPILE_SMP
    valid_config = valid_config && config.contains("num_processors");
#endif
    
    if (valid_config) {
      // Read all mandatory config file settings.
      inseq1         = config.getOption<string>("inseq1");
      inseq2         = config.getOption<string>("inseq2");
      outct          = config.getOption<string>("outct");
      outct2         = config.getOption<string>("outct2");
      aout           = config.getOption<string>("aout");
      
	  

#ifdef COMPILE_SMP
      numProcessors = config.getOption<int>("num_processors");
#endif
    }

    // Read optional config file settings
    if (valid_config) {

		//check to see if a ct file is specified for sequence 1, i.e. ct_templated = true.
	   if (config.contains("ct_templated")) ct_templated = config.getOption<bool>("ct_templated");
	   else ct_templated = false;
	  if (config.contains("imaxseparation")) imaxseparation = config.getOption<int>("imaxseparation");
	  else imaxseparation = -99;//set a flag to indicate the M is not being used

#ifdef DYNALIGN_II
          if (config.contains("slope")) slope           = config.getOption<double>("slope");
	  else slope = 0.1;
          if (config.contains("intercept")) intercept           = config.getOption<double>("intercept");
	  else intercept = 0.5;
          if (config.contains("max_elongation")) max_elongation           = config.getOption<int>("max_elongation");
	  else max_elongation = 5;
#else
#endif
          if (config.contains("fgap")) fgap           = config.getOption<double>("fgap");
	  else fgap = 0.4;
      if (config.contains("maxtrace")) maxtrace       = config.getOption<int>("maxtrace");
	  else maxtrace = 750;
      if (config.contains("percent")) percent        = config.getOption<int>("percent");
	  else percent =20;
      if (config.contains("singlefold_subopt_percent")) singlefold_subopt_percent =
                       config.getOption<int>("singlefold_subopt_percent");
	  else singlefold_subopt_percent=30;
      if (config.contains("bpwin")) bpwin          = config.getOption<int>("bpwin");
	  else bpwin = 2;
      if (config.contains("awin")) awin           = config.getOption<int>("awin");
	  else awin = 1;
#ifdef DYNALIGN_II
#else
      if (config.contains("insert")) insert         = config.getOption<bool>("insert");
	  else insert = true;
#endif
      if (config.contains("dsv_templated")) dsv_templated  = config.getOption<bool>("dsv_templated");
	  else dsv_templated = false;
      if (config.contains("savefile")) {
        savefile = config.getOption<string>("savefile");
      }
      if (config.contains("maxpairs")) maxpairs = config.getOption<int>("maxpairs");
		else maxpairs = -1;
      if (config.contains("constraint_1_file")) {
        constraint1 = config.getOption<string>("constraint_1_file");
      }

      if (config.contains("constraint_2_file")) {
        constraint2 = config.getOption<string>("constraint_2_file");
      }

      if (config.contains("shape_1_file")) {
        shape1 = config.getOption<string>("shape_1_file");
      }

	  if (config.contains("shape_2_file")) {
        shape2 = config.getOption<string>("shape_2_file");
      }

	  if (config.contains("constraint_align_file")) {
        constrainta = config.getOption<string>("constraint_align_file");
      }

	  if (config.contains("optimal_only")) {
		optimalonly = config.getOption<bool>("optimal_only");
	  }
	  //by default, setup to determine suboptimal structures
	  else optimalonly=false;
	  if (config.contains("local")) {
		local = config.getOption<bool>("local");
	  }
	  //by default, run global alignment
	  else local=false;

	  //Read SHAPE parameters or set defaults
	  if (config.contains("shapeslope1")) shapeslope1 = config.getOption<double>("shapeslope1");
	  else shapeslope1 = 1.8;
	  if (config.contains("shapeslope2")) shapeslope2 = config.getOption<double>("shapeslope2");
	  else shapeslope2 = 1.8;
	  if (config.contains("shapeintercept1")) shapeintercept1 = config.getOption<double>("shapeintercept1");
	  else shapeintercept1 = -0.6;
	  if (config.contains("shapeintercept2")) shapeintercept2 = config.getOption<double>("shapeintercept2");
	  else shapeintercept2 = -0.6;

	  if (config.contains("maximumpairingdistance1")) maximumpairingdistance1 = config.getOption<int>("maximumpairingdistance1");
	  else maximumpairingdistance1 = 0;
	  if (config.contains("maximumpairingdistance2")) maximumpairingdistance2 = config.getOption<int>("maximumpairingdistance2");
	  else maximumpairingdistance2 = 0;
	  if (config.contains("temperature")) Temperature = config.getOption<double>("temperature");
	  else Temperature = 310.15;
	  if (config.contains("DNA")) DNA = config.getOption<bool>("DNA");
	  else DNA = false;
	  // select the extended alphabet or use "rna" or "dna" depending on the DNA bool value.
	  if (config.contains("alphabet")) alphabet = config.getOption<string>("alphabet");
	  else alphabet = DNA?DT_DNA:DT_RNA;
    }

    // Read settings dependent on previous settings.
    if (valid_config && dsv_templated) {
      valid_config =
        config.contains("dsvtemplatename") &&
        config.contains("maxdsvchange");
      if (valid_config) {
        dsvtemplatename = config.getOption<string>("dsvtemplatename");
        maxdsvchange    = config.getOption<float>("maxdsvchange");
      }
    }
    
    if (!valid_config) {
      cerr << "ERROR: At least one parameter could not be read from the "
           << "configuration file. Aborting." << endl;
      exit(1);
    }
      
  } else {
	  
	  //set ct_templated to false by default:
	  ct_templated = false;
      
    // Interactive mode:
#ifdef DYNALIGN_II
	cout << "Usage: dynalign_ii [config file]\nNo config file specified, using interactive mode:\nNote that most, but not all, features are available through the interactive mode.\n\n";
#else
	cout << "Usage: dynalign [config file]\nNo config file specified, using interactive mode:\nNote that most, but not all, features are available through the interactive mode.\n\n";
#endif      
    cout << "Enter the name of the first sequence: ";
    cin >> inseq1;
	    
    cout << "Enter the name of the second sequence: ";
    cin >> inseq2;
	    
    cout << "Enter the name of the first output ct file: ";
    cin >> outct;
	    
    cout << "Enter the name of the second output ct file: ";
    cin >> outct2;
	    
    cout << "Enter name for the output of the alignment: ";
    cin >> aout;
	    
    cout << "Enter the max separation (M) <Recommend enter -99 for constraint by HMM forward-backward>: ";
    cin >> imaxseparation;
#ifdef DYNALIGN_II
    cout << "Enter the slope of domain insertion penalty: <Recommend 0.1>";
    cin >> slope;

    cout << "Enter the intercept of domain insertion penalty: <Recommend 0.5>";
    cin >> intercept;


    cout << "Enter the maximum value for helix elongation: <Recommend 5>";
    cin >> max_elongation;
#else
#endif
    cout << "Enter the gap penalty: <Recommend 0.4>";
    cin >> fgap;
	    
    cout << "Enter the maximum number of structures: ";
    cin >> maxtrace;
	    
    cout << "Enter the maximum percent energy difference "
         << "(where entering 20 is 20%): ";
    cin >> percent;

    cout << "Enter the maximum percent energy difference for keeping "
         << "suboptimal structures from the single-sequence folding "
         << "calculations: <Recommend 30>";
    cin >> singlefold_subopt_percent;
	    
    cout << "Enter the base pair window size: <Recommend 2>";
    cin >> bpwin;
	    
    cout << "Enter the alignment window size: <Recommend 1>";
    cin >> awin;
#ifdef DYNALIGN_II
#else    
    cout << "Allow single BP inserts into one sequence? (1/0) <Recommend 1>";
    cin >> insert;
#endif	    
    cout << "Write save file (needed for refold or dot plots)? (1/0) ";
    cin >> query;
	    
    if (query) {
      cout << "Enter the save file name: ";
      cin >> savefile;
    }
	    
    cout << "Input constraints for sequence 1? (1/0) ";
    cin >> query;
	    
    if (query) {
      cout << "Enter the constraint file name: ";
      cin >> constraint1;
    }
	    
    cout << "Input constraints for sequence 2? (1/0) ";
    cin >> query;
	    
    if (query) {
      cout << "Enter the constraint file name: ";
      cin >> constraint2;
    }
	    
    cout << "Input constraints for alignment? (1/0) ";
    cin >> query;
	    
    if (query) {
      cout << "Enter the constraint file name: ";
      cin >> constrainta;
    }
    
    cout << "Dsv template to be used (i.e. a progressive calculation) "
         << "(1/0): ";
    cin >> dsv_templated;
    if(dsv_templated)
      {
        cout << "Specify dsv template file: ";
        cin >> dsvtemplatename;
        cout << "Percent energy window for bases to keep from dsv template "
             << "(where 20 = 20% difference from optimal): ";
        cin >> maxdsvchange;
		cout << "Number of base pairs to allow for template <recommend -1, meaning the length of the 1st sequence: ";
		cin >> maxpairs;
      }
#ifdef COMPILE_SMP
    cout << "Enter the number of processors on this machine: ";
    cin >> numProcessors;
#endif
  }

  // Validate settings
#ifdef DYNALIGN_II
  islope = (int)(slope * 10.0);
  iintercept = (int)(intercept * 10.0);
#else
#endif

  igapincrease = (int)(fgap * 10.0);

#ifdef COMPILE_SMP
  if (numProcessors < 1) {
    cerr << "Warning: invalid number of processors specified; "
         << "defaulting to 1." << endl;
    numProcessors = 1;
  }
#endif

//cdbg << "After config" << endl;

	//open the thermodynamic data tables
  Thermodynamics thermo(!DNA, alphabet.c_str(), Temperature);
  thermo.ReadThermodynamic();
	RNA rna1(inseq1.c_str(), ct_templated?FILE_CT:FILE_SEQ, &thermo); // ct_templated means the first sequence should be read from a ct file, not a .seq file.
	RNA rna2(inseq2.c_str(), FILE_SEQ, &thermo);
	ErrorChecker<RNA> ec1(&rna1), ec2(&rna2);
	if (ec1.isErrorStatus()||ec2.isErrorStatus())
		return 1;

  //cdbg << "After RNA" << endl;

	structure &ct1 = *rna1.GetStructure(), &ct2 = *rna2.GetStructure();

#ifdef CHECK_ARRAY
        seq1 = inseq1;
        seq2 = inseq2;
#endif
	//allocate space for the alignment
	align = new short *[maxtrace];//maximum number of tracebacks and next line and below at delete
	for (i=0;i<maxtrace;i++)  align[i] = new short [rna1.GetSequenceLength()+1];

	//do the dynalign calculation
  //cdbg << "Before constraints" << endl;
	constrained = false;
	if (constraint1 != "") {
		if (ec1.isErrorStatus(rna1.ReadConstraints(constraint1.c_str()))) return 1;
		constrained = true;
	}
	if (constraint2 != "") {
		if (ec2.isErrorStatus(rna2.ReadConstraints(constraint2.c_str()))) return 1;
		constrained = true;
	}
	
	if (shape1 != "") {
		if (ec1.isErrorStatus(rna1.ReadSHAPE(shape1.c_str(), shapeslope1, shapeintercept1))) return 1;
		constrained = true;		
	}
	if (shape2 != "") {
		if (ec2.isErrorStatus(rna2.ReadSHAPE(shape2.c_str(), shapeslope2, shapeintercept2))) return 1;
		constrained = true;		
  }
  
  //cdbg << "Before constraints Aln" << endl;
	if (constrainta != "") {
		//Apply alignment constraints

		// check that the file exists.
		if (!fileExists(constrainta, true)) {
			cerr << "Alignment Constraints file, \""<<constrainta<<"\" does not exist or could not be opened." << endl;
			return 1;
		}
		constrained = true;
		forcealign=new short *[2];
		
		forcealign[0]=new short [rna1.GetSequenceLength()+1];
		forcealign[1]=new short [rna2.GetSequenceLength()+1];
		for (index=1;index<=rna1.GetSequenceLength();index++) {
			forcealign[0][index]=0;
		}
		for (index=1;index<=rna2.GetSequenceLength();index++) {
			forcealign[1][index]=0;
		}
		readalignmentconstraints(constrainta.c_str(),forcealign,&ct1,&ct2);
	} 
	else {
		forcealign = NULL;
	}
  
	//check if lowercase nucleotides were entered in either sequences - these will be single-stranded
	if (ct1.GetNumberofSingles()>0||ct2.GetNumberofSingles()>0) constrained=true;
  
  //cdbg << "After constraints" << endl;
	
	
	//This section folds thie individual sequences to find insignificant pairs to be ingnored by Dynalign.
	ct1.allocatetem();
	ct2.allocatetem();

	datatable &data = *thermo.GetDatatable();

	//dsv_templated is helpful for progressive calculations from multiple dynalign calculations
	if(dsv_templated) templatefromdsv( &ct1, dsvtemplatename.c_str(), maxdsvchange, maxpairs);

	//ct_templted is helpful if the structure is known for sequence 1
	else if (ct_templated) templatefromct(&ct1);

	//otherwise, fold the sequence to determine which pairs can result in low free energy structures
	else { 

		//if maximumpairingdistance1 is > 0 (the default), apply this to the templatefromfold calulation
		if (maximumpairingdistance1>0) {
			ct1.SetPairingDistance(maximumpairingdistance1);
		}
		templatefromfold(&ct1, &data, singlefold_subopt_percent);

	}

	//if maximumpairingdistance2 is > 0 (the default), apply this to the templatefromfold calulation
	if (maximumpairingdistance2>0) {
		
		ct2.SetPairingDistance(maximumpairingdistance2);
	}
	templatefromfold( &ct2, &data, singlefold_subopt_percent );

	//This next section determined the allowed nucleotide alignments if the HMM forward-backward is used:
	if (imaxseparation == -99) {
		//allocate space in allowed_alignments
		allowed_alignments = new bool *[ct1.GetSequenceLength()+1];
		for (i=0;i<=ct1.GetSequenceLength();i++) {
			allowed_alignments[i] = new bool [ct2.GetSequenceLength()+1];	
	
		}

		//calculate_aln_probs_env(&ct1, &ct2, NULL, allowed_alignments, align_threshold);
		calculate_coinc_probs_env(&ct1, &ct2, allowed_alignments, forcealign);
	}
	else allowed_alignments = NULL;

	int resultcode;
#ifdef DYNALIGN_II

  resultcode = dynalign(&ct1, &ct2, align, imaxseparation, islope, iintercept, igapincrease, &data,
                 maxtrace, bpwin, awin, percent, forcealign, max_elongation, allowed_alignments, &progress,
           (savefile != "") ? savefile.c_str() : NULL, optimalonly, local,
		   /*force =*/ constrained, numProcessors);
#else
  resultcode = dynalign(&ct1, &ct2, align, imaxseparation, igapincrease, &data,
           insert, maxtrace, bpwin, awin, percent, forcealign, allowed_alignments, &progress,
           (savefile != "") ? savefile.c_str() : NULL, optimalonly, local,
		   /*force =*/ constrained, numProcessors);
#endif

  if (resultcode==14)
		//dynalign returns an int that indicates an error.  14 is currently the only possible error, a traceback error
		cerr << "Dynalign encountered a traceback error.  Please report this possible bug to David_Mathews@urmc.rochester.edu" << endl;
  else if (resultcode!=0)
		cerr << "Dynalign encountered an unknown error #" << resultcode << ".  Please report this possible bug to David_Mathews@urmc.rochester.edu" << endl;

	if (resultcode!=0) return 1;
	
	//output the structures
	int err1 = ec1.isErrorStatus(rna1.WriteCt(outct.c_str()));
	int err2 = ec2.isErrorStatus(rna2.WriteCt(outct2.c_str()));

	//output the alignment
	alignout(align,aout.c_str(),&ct1,&ct2);

	if (err1 || err2) return 1;

	//clean up memory allocation
	for (i=0;i<maxtrace;i++)  delete[] align[i];
	delete[] align;

	if (constrainta != "") {
		delete[] forcealign[0];
		delete[] forcealign[1];
		delete[] forcealign;
	}

	if (imaxseparation == -99) {
		//delete space in allowed_alignments
		for (i=0;i<=ct1.GetSequenceLength();i++) {
			delete[] allowed_alignments[i];	
		}
		delete[] allowed_alignments;
	}
	return 0;
}


