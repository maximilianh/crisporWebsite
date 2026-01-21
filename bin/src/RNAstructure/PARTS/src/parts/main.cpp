#include "ppf_cli.h"
#include "ppf_loops.h"

#include <iostream>
using namespace std;

int main( int argc, char* argv[] ) {

  // Do error checking on the command line, make sure its usage is correct.
  if( argc != 2 ) {
    cerr << "Incorrect number of command line arguments given." << endl
         << "USAGE: PARTS <configuration file>" << endl << endl
	 << "Configuration file format:" << endl
	 << "<flag> <value> <new line> ..." << endl << endl
	 << "Configuration file flags:" << endl
	 << "seq1" << endl
	 << "    First input sequence file."
	 << endl << endl
	 << "seq2" << endl
	 << "    Second input sequence file."
	 << endl << endl
	 << "mode" << endl
	 << "    Mode used during calculation." << endl
	 << "    Value can be map, pp, or stochsample (stochastic sampling)."
	 << endl << endl
	 << "seq1_map_ct_op" << endl
	 << "    Output CT file name for sequence 1, map mode."
	 << endl << endl
	 << "seq2_map_ct_op" << endl
	 << "    Output CT file name for sequence 2, map mode."
	 << endl << endl
	 << "map_aln_op" << endl
	 << "    Alignment output file name, map mode."
	 << endl << endl
	 << "seq1_pp_op" << endl
	 << "    Probability matrix output file name for sequence 1, pp mode."
	 << endl << endl
	 << "seq2_pp_op" << endl
	 << "    Probability matrix output file name for sequence 2, pp mode."
	 << endl << endl
	 << "seq1_sample_ct_op" << endl
	 << "    Output CT file for sequence 1, stochastic sampling mode."
	 << endl << endl
	 << "seq2_sample_ct_op" << endl
	 << "    Output CT file for sequence 2, stochastic sampling mode."
	 << endl << endl
	 << "sample_aln_op" << endl
	 << "    Alignment output file name, stochastic sampling mode."
	 << endl << endl
	 << "nsamp" << endl
	 << "    Number of sampled structures to generate." << endl
	 << "    (stochastic sampling mode only)"
	 << endl << endl
         << "seed" << endl
         << "    Pseudo-random number generator seed." << endl
         << "    (stochastic sampling mode only)"
         << endl << endl
	 << "Flags 'seq1,' 'seq2', and 'mode' are required to run PARTS."
	 << endl
	 << "The other flags are optional." << endl << endl
	 << "Note that if output file names are not specified, output files "
	 << "are instead"
	 << endl
	 << "placed in the current directory, with default names."
	 << endl;
    return 0;
  }

  // Create a variable that holds the error state, and another variable that
  // holds whether a PARTS object was created successfully.
  bool isError = false;
  bool partsCreated = false;

  // Initialize a t_ppf_loops object with the command line arguments.
  // Then, check to make sure the object was created successfully.
  // If it was created successfully, continue. If not, show an error.
  cout << "Initializing PARTS...." << flush;

  t_ppf_loops* PARTS = new t_ppf_loops( argc, argv );
  if( PARTS->GetErrorCode() != 0 ) {
    cerr << endl << PARTS->GetErrorMessage( PARTS->GetErrorCode() ) << endl;
    isError = true;
  }

  if( !isError ) { cout << "done." << endl; }

  // Run the PARTS calculation with the t_ppf_loops object. Only do this
  // calculation if the t_ppf_loops object was created successfully.
  // Then, check to make sure the calculation ran successfully.
  // If it ran successfully, continue. If not, show an error.
  if( PARTS->GetErrorCode() == 0 ) {
    cout << "Computing main calculation..." << endl;

    PARTS->compute_pseudo_free_energy_loops();
    if( PARTS->GetErrorCode() != 0 ) {
      cerr << endl << PARTS->GetErrorMessage( PARTS->GetErrorCode() ) << endl;
      isError = true;
    } else { partsCreated = true; }

    if( !isError ) { cout << endl << "done." << endl; }
  }

  // Delete PARTS object, if necessary.
  if( partsCreated ) { delete PARTS; }

  // Print out a confirmation of run finishing and return 0.
  if( !isError ) { cout << "PARTS complete." << endl; }
  else { cerr << "PARTS complete with errors." << endl; }

  return 0;
}

