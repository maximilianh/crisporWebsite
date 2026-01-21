//
// Created by Max Ward on 27 March 2017.
//


/**
 * This file defines a class that implements that nearest neighbour energy model as described by Mathews & Turner,
 * Zuker et al., and Tinoco et al.
 */

#ifndef NN_MODEL_HPP
#define NN_MODEL_HPP

#include "primary_structure.hpp"
#include "energy.hpp"

#include "../src/rna_library.h"
#include "../src/structure.h"

/**
 * An instance of NNModel represents a particular parametrization of the nearest neighbour energy model.
 * Currently based largely on the 2004 version found at http://rna.urmc.rochester.edu/NNDB/.
 * Because multi-loops are treated oddly (linear for DP, logarithmic for energy evaluation) they are abstracted out.
 * Extensions of this class are used to provide energy parameters and functions for different multi-loop energy models.
 * Note that all the energy functions are orthogonal, so, for example, no function will call Branch for you.
 * In addition, the functions provided have a philosophy of being as clever as possible. They will deduce what kind of
 * coaxial stack you mean--whether i,j closes a multi-loop--for example.
 *
 * Note that many of these functions make documented assumptions. These are generally checked by standard C assertions,
 * which can be turned off using the NDEBUG flag.
 *
 * Note that zero indexing is assumed for base indexes!
 */
class NNModel {
protected:
	efnmax::PrimeStructure rna;
	datatable *dt = NULL;
	structure *struc = NULL;
public:
	/// Minimum number of unpaired nucleotides allowed in a hairpin loop.
	static const int MIN_HAIRPIN_UNPAIRED = 3;

	/**
	 * Returns the energy of a one-loop from 5' i to 3' j.
	 * A one-loop is a hairpin loop.
	 * Assumes i<j.
	 * http://rna.urmc.rochester.edu/NNDB/turner04/hairpin.html
	 * @param i The 5' nucleotide.
	 * @param j the 3' nucleotide.
	 * @return The free energy change.
	 */
	virtual efnmax::energy_t OneLoop(int i, int j) const;


	/**
	 * Returns the energy of a two-loop between 5' -> i -> k -> l -> j -> 3'.
	 * (_(_)_)
	 * i k l j
	 * This could be a bulge, stacking region, or internal loop.
	 * http://rna.urmc.rochester.edu/NNDB/turner04/bulge.html
	 * http://rna.urmc.rochester.edu/NNDB/turner04/wc.html (stacking)
	 * http://rna.urmc.rochester.edu/NNDB/turner04/internal.html
	 * @param i 5' End of the closing bond.
	 * @param k 5' End of the internal bond.
	 * @param l 3' End of the internal bond.
	 * @param j 3' End of the closing bond.
	 * @return The free energy change. Possibly a large 'infinity' value if the configuration is invalid.
	 */
	virtual efnmax::energy_t TwoLoop(int i, int k, int l, int j) const;


	/**
	 * The free energy cost of a (external/multi)-loop branch closed by 5' i and 3' j.
	 * Considers AU/GU closure penalties.
	 * Note that this function doesn't care about i<j.
	 * @param i 5' end of the branch.
	 * @param j 3' end of the branch.
	 * @return The free energy change.
	 */
	virtual efnmax::energy_t Branch(int i, int j) const;


	/**
	 * A flush coaxial stacking between pairs i,j and k,l.
	 * 5' i and 3' j are a bonding pair which are in a coaxial stack with 5' k and 3' l.
	 * Assumes i<j and k<l.
	 * Assumes i,j < k,l, or that i,j closes a multi-loop. The function will deduce this information.
	 * Assumes that i,j and k,l are directly adjacent branches.
	 * Note that this only computes the coaxial stacking score.
	 * Thus, Branch WILL NOT be called by this function.
	 * http://rna.urmc.rochester.edu/NNDB/turner04/coax.html
	 * @param i 5' end of the first branch.
	 * @param j 3' end of the first branch.
	 * @param k 5' end of the second branch.
	 * @param l 3' end of the second branch.
	 * @return The free energy change of the stack. Possibly a large 'infinity' value if the configuration is invalid.
	 */
	virtual efnmax::energy_t FlushCoax(int i, int j, int k, int l) const;


	/**
	 * Mismatch mediated coaxial stack. Mismatched bases are next to stem at i,j.
	 * Parameters k and l represent the other stem in the coaxial stack.
	 * Assumes i<j and k<l.
	 * Assumes i,j < k,l, or that i,j or k,l close a multi-loop. The function will deduce this information.
	 * Assumes i,j and k,l are branches separated by exactly one nucleotide.
	 * Note that this only computes the coaxial stacking score.
	 * Thus, Branch WILL NOT be called by this function.
	 * http://rna.urmc.rochester.edu/NNDB/turner04/coax.html
	 * @param i 5' end of the first branch.
	 * @param j 3' end of the first branch.
	 * @param k 5' end of the second branch.
	 * @param l 3' end of the second branch.
	 * @return The free energy change of the stack. Possibly a large 'infinity' value if the configuration is invalid.
	 */
	virtual efnmax::energy_t MismatchCoax(int i, int j, int k, int l) const;


	/**
	 * Returns the free energy change of a 5' dangle off the pair i,j.
	 * Assumes i<j and that i,j DO NOT close a multi-loop. Which means dangle<i
	 * @param i 5' end of the branch.
	 * @param j 3' end of the branch.
	 * @return The free energy change. Possibly a large 'infinity' value if the configuration is invalid.
	 */
	virtual efnmax::energy_t FiveDangle(int i, int j) const;

	/**
	 * Returns the free energy change of a 5' dangle off the pair i,j.
	 * Assumes i<j and that i,j DO close a multi-loop. Which means i<dangle<j.
	 * @param i 5' end of the branch.
	 * @param j 3' end of the branch.
	 * @return The free energy change. Possibly a large 'infinity' value if the configuration is invalid.
	 */
	virtual efnmax::energy_t ClosingFiveDangle(int i, int j) const;

	/**
	 * Returns the free energy change of a 3' dangle off the pair i,j.
	 * Assumes i<j and that i,j DO NOT close a multi-loop. Which means j<dangle.
	 * @param i 5' end of the branch.
	 * @param j 3' end of the branch.
	 * @return The free energy change. Possibly a large 'infinity' value if the configuration is invalid.
	 */
	virtual efnmax::energy_t ThreeDangle(int i, int j) const;

	/**
	* Returns the free energy change of a 3' dangle off the pair i,j.
	* Assumes i<j and that i,j DO close a multi-loop. Which means i<dangle<j.
	* @param i 5' end of the branch.
	* @param j 3' end of the branch.
	* @return The free energy change. Possibly a large 'infinity' value if the configuration is invalid.
	*/
	virtual efnmax::energy_t ClosingThreeDangle(int i, int j) const;

	/**
	 * Returns the free energy change of a terminal mismatch off the pair i,j.
	 * Assumes i<j and that i,j DO NOT close a multi-loop. Which means mismatch5'<i<j<mismatch3'.
	 * @param i 5' end of the branch.
	 * @param j 3' end of the branch.
	 * @return The free energy change. Possibly a large 'infinity' value if the configuration is invalid.
	 */
	virtual efnmax::energy_t Mismatch(int i, int j) const;

	/**
	 * Returns the free energy change of a terminal mismatch off the pair i,j.
	 * Assumes i<j and that i,j DO close a multi-loop. Which means i<3'dangle<5'dangle<j.
	 * @param i 5' end of the branch.
	 * @param j 3' end of the branch.
	 * @return The free energy change. Possibly a large 'infinity' value if the configuration is invalid.
	 */
	virtual efnmax::energy_t ClosingMismatch(int i, int j) const;

	/// Returns a safe "infinity" value.
	virtual efnmax::energy_t MaxMFE() const;


	/**
	 * Sets the current working RNA for the energy functions.
	 * @param rna The new working RNA.
	 */
	void SetRNA(const std::string &rna);


	/// Assumes the caller will manage the data table memory.
	NNModel(datatable *_dt, const std::string &_rna) {
		this->dt = _dt;
		SetRNA(_rna);
	}

	virtual ~NNModel() {
		delete this->struc;
	}

};

#endif //EFNMAX_NN_MODEL_HPP
