//
// Created by Max Ward on 29 March 2017.
//

#ifndef EFNMAX_NN_AFFINE_MODEL_HPP
#define EFNMAX_NN_AFFINE_MODEL_HPP

#include "nn_model.hpp"
#include "energy.hpp"

/**
 * Extends the nearest neighbour model to include a linear energy function for multi-branch loops.
 * This is the model typically employed by dynamic programming algorithms for RNA folding.
 */
class NNModelLinear: public NNModel {

public:
	/**
	 * @return The energy of multi-loop closure under the linear model.
	 * @param branches The number of branches in the multi-loop. Assumed to be >= 3.
	 * @param unpaired The number of unpaired nucleotides in the multi-loop.
	 */
	efnmax::energy_t MLClosure(unsigned branches, unsigned unpaired) const;
	/**
	 * The initiation cost of a multi-loop.
	 */
	efnmax::energy_t MLInitCost() const;
	/**
	 * The bonus given the a branch in a multi-loop.
	 */
	efnmax::energy_t MLBranchCost() const;
	/**
	 * The penalty given to an unpaired nucleotide in a mutli-loop.
	 */
	efnmax::energy_t MLUnpairedCost() const;
	NNModelLinear(datatable *_dt, const std::string &_rna)
		: NNModel(_dt, _rna) {}
};

#endif //EFNMAX_NN_AFFINE_MODEL_HPP
