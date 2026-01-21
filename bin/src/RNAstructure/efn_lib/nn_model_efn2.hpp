//
// Created by Max Ward on 29 March 2017.
//

#ifndef EFNMAX_NN_EFN2_MODEL_HPP
#define EFNMAX_NN_EFN2_MODEL_HPP

#include "nn_model_linear.hpp"
#include "energy.hpp"

/**
 * Extends the nearest neighbour model using the linear model to include the energy function used by EFN2 to score multi-loops.
 * This model is the same as for EFN2 from RNAstructure 5.8.
 * The model has the form: 
 * mlasym*max(average_asymmetry, 2) + efn2a + efn2c*#branches + efn2b*#unpaired_nucleotides + strain if #unpaired_nucleotides <= 8
 * mlasym*max(average_asymmetry, 2) + efn2a + efn2c*#branches + efn2b*8 + 1.1*ln(#unpaired_nucleotides/8) + strain if #unpaired_nucleotides > 8
 * See this (https://rna.urmc.rochester.edu/NNDB/turner99/mb.html) for information about terms not relating to asymmetry or strain.
 * See this (https://rna.urmc.rochester.edu/NNDB/turner04/mb.html) for the asymmetry or strain terms.
 */
class NNModelEFN2: public NNModelLinear {
private:
	int ml_pivot;
	efnmax::kcalmol_t ml_log_multiplier;

public:
	/**
	 * @return The energy of multi-loop closure under the EFN2 model.
	 * @param average_asymmetry The average asymmetry of unpaired nucleotides on either side of each branch.
	 * @param branches The number of branches in the multi-loop. Assumed to be >= 3.
	 * @param unpaired The number of unpaired nucleotides in the multi-loop.
	 */
	efnmax::energy_t MLClosure(double average_asymmetry, unsigned unpaired, unsigned branches) const;
	/**
	 * The initiation cost of a multi-loop.
	 */
	efnmax::energy_t MLAverageAsymmetryCost() const;
	/**
	 * The penalty for a strained multi-branch.
	 */
	efnmax::energy_t MLStrainCost() const;

	NNModelEFN2(datatable *_dt, const std::string &_rna)
		: NNModelLinear(_dt, _rna) {
			ml_pivot = 8;
			ml_log_multiplier = 1.1;
		}
};

#endif //EFNMAX_NN_EFN2_MODEL_HPP
