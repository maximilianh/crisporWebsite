#if !defined(INTERMOLECULAR_H)
#define INTERMOLECULAR_H

/*
intermolecular.h and intermolecular.cpp inculde funcions calculating and report
different free energy for binding in OligoWalk.
They are revised based on Mathews' code from RNAStructure.
olig() generate the energy data;
report save the generated data;
siprefileter and filterbysirna() were used to prefilter and postfileter the 
functional siRNA respectively, using criterias other than thermodynamics


															----Nov., 2005
															John(Zhi Lu)
# Updates to report output 2018-01-18 Richard Watson.
*/

#include <cmath>
//#include "stdafx.h"
#include "rna_library.h"
#include "thermo.h"
#include "pclass.h"
#include "algorithm.h"
#include "siRNAfilter.h"
#include "alltrace.h"
#include "alltrace_intermolecular.h"
#include "stochastic.h"
#include "TProgressDialog.h"
#include "../src/phmm/utils/xmath/log/xlog_math.h"


#define FILTER_PASS 6	//The criteria for prefiltering of functional siRNA in OligoWalk
#define maxfile 75		//maximum length of file names

//this structure rddata contains all the thermodynamic parameters needed
//for determining the stability of an RNA-DNA duplex
struct rddata {

   short int stack[5][5][5][5];
   short int init;
};


//! olig is the backend function for oligo walk

//! table	-- A table (2D array) that is filled with the thermodynamic data 
//! 	for each oligo. The first dimension is an array of rows, 
//!		each representing an oligo bound to a specific site on the target.
//! 	The second dimension is an array with the following indices:
//! 	[0] - overall DG
//! 	[1] - duplex DG
//! 	[2] - free energy of breaking target structure
//! 	[3] - intramolecular oligo free energy
//! 	[4] - intermolecular oligo free energy
//! 	[5] - Tm of duplex (x10)
//! 
//! mode	-- How the target structure should be handled.
//! 	1 - break local target structure to bind oligo 
//! 	2 - refold target RNA after oligo binding
//! 	3 - no target structure considered
//! 
//! suboptimal - How suboptimal structures should be handled.
//! 	0 - only consider optimal structure
//! 	1 - like choice 3,using suboptimal structures,but the whole set from alltarce() function prediction
//! 	2 - using partition funcion considering every possible structure  of target
//! 		suboptimal 2 can only used with mold 2
//! 	3 - using suboptimal structures (heuristic method) for both oligo-free and oligo-bound target RNA
//! 	4 - using stochasit sampling method to sample 1000 structures
//! 
//! prefilter 	-- Whether criteria should be used to prefill functional siRNA
//! 	0 - no prefilter
//! 	1 - use prefilter
//! 
//! foldsize	-- Only fold a fragment with size=foldsize+binding_length, 
//! 	which is centered on the siRNA binding region.
//! 	when foldsize > 1, only mode=2 and suboptimal=2 are valid options.
//! 
//! distance 	-- limit the maximum distance between nucleotides that can pair
//! 
//! shapefile	-- specify a SHAPE datafile (set to "" to not use SHAPE data)
//! 
//! TESTS	-- run tests, set to -1 for no tests.
//!
//! writesav - write sav files to save time in test mode
//!
void olig(const bool isdna, const int mode, structure *ct, const int length, const double conc, const int suboptimal,
	const int start, const int stop, int foldsize, const int distance, 
	int **table, int **numofsubstructures, const char *const shapefile, int *TEST, const bool writesav,
	datatable& data, datatable& ddata, thermo* helixstack, rddata *hybriddata, 
	siPREFILTER *prefilter, ProgressHandler *update);


//this function reads data into structure rddata
int readrd (rddata* data, const string &dnarnafile);

//! Function: report(ostream& out, ...)
//! This function writes a tab delimited (or HTML) file with the oligowalk data
//! Parameters:
//! out			-- The output stream (e.g. to an output file).
//! ct			-- Contains the sequence
//! table		-- The table filled by OligoWalk
//! numofsubstructure	-- 
//! length		-- The length of the oligonucleotides
//! isdna		-- Whether the oligos are dna (true = DNA, false = RNA)
//! conc		-- The oligonucleotide concentration
//! suboptimal	-- Whether suboptimal structures were used
//! start		-- The nucleotide postition of the start of the walk on the target sequence
//! stop		-- The nucleotide postition of the end of the walk on the target sequence
//! prefilter	-- John Lu's siRNA prefilter
//! foldsize	-- The size of the folding region that is centered around the oligo
//! mask		-- A bool array in which mask[i] is true if the oligo at table[i] has passed 
//! 			   the filter criteria (i.e. asuf, tofe, fnnfe).
//! asuf		-- Minimum Antisense strand unimolecular folding free energy
//! tofe		-- Minimum Target opening free energy
//! fnnfe		-- Minimum First nearest neighbor free energy
//! isHTML		-- Whether this is HTML or tab delimited text (false = tab delimited, true = HTML)
//! writeHeader	-- Whether header information should be written in the report.
//! writeBody	-- Whether the report body should be written (as opposed to just the header).
void report(ostream& out, 
	const bool isdna, const int mode, structure *ct, const int length, const double conc, const int suboptimal,
	const int start, const int stop, const int foldsize, const int distance,
	int **table, int **numofsubstructures, const char * const shapefile,
	const siPREFILTER *prefilter,
	const bool *mask=NULL, const double asuf=0, const double tofe=0, const double fnnfe=0,
	const bool isHTML=true, const bool writeHeader=true, const bool writeBody=true);

// void report(ostream& out, structure *ct, int **table,int **numofsubstructures, const int length, const bool isdna,
// 			const double conc, const int suboptimal,const int start,const int stop,const siPREFILTER *prefilter,const int foldsize,
// 			const bool *mask=NULL, const double asuf=0, const double tofe=0, const double fnnfe=0, const bool isHTML=true);

//returns the base complementary to base i
//with 1 == A, 2 == C, 3 == G, 4 == U
int complement(const int i, structure *ct);

//returns the character (ACGT/U) corresponding to the base-number.
char numtobase(const int basenum, structure *ct, const bool isDNA);

//post-filter of functional siRNA written by Dave.
void filterbysirna(structure *ct, int **table, int length, datatable *data, 
				   bool *mask, double asuf, double tofe, double fnnfe);

//copy the arrays to be reused with different index when the folded region move one nucleotide to the right
inline void scancopy(OligoPclass *region, OligoPclass *copyregion);
//copy every arrays to be reused when the folded region did not move right yet 
inline void scancopyend(OligoPclass *region, OligoPclass *copyregion) ;

#endif

