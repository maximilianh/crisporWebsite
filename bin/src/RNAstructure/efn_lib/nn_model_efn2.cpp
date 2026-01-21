//
// Created by Max Ward on 31 March 2017.
//

#include "nn_model_efn2.hpp"

#include <cmath>
#include <algorithm>
#include <cassert>

#include <iostream>

efnmax::energy_t NNModelEFN2::MLClosure(double average_asymmetry, unsigned branches, unsigned unpaired) const {
	assert(branches >= 3);

	efnmax::energy_t branch_cost = branches * MLBranchCost();
	efnmax::kcalmol_t logarithmic_part;
	if (unpaired <= ml_pivot)
		logarithmic_part = efnmax::EnergyToKCal(MLInitCost() + branch_cost + unpaired * MLUnpairedCost());
	else // Use D*(int)round(log(...)) to get same behaviour as RNAstructure.
		logarithmic_part = efnmax::EnergyToKCal(MLInitCost() + branch_cost + ml_pivot * MLUnpairedCost())
			+ ml_log_multiplier * std::log(unpaired / (double) ml_pivot);

	efnmax::kcalmol_t strain_part = efnmax::EnergyToKCal(branches == 3 && unpaired < 2 ? MLStrainCost() : 0);

	// Needs to convert to kcalmol_t manually to ensure no precision is lost, even though MLAverageAsymmetryCost is energy_t.
	efnmax::kcalmol_t asymmetry_part = std::min(2.0, average_asymmetry)*MLAverageAsymmetryCost()*efnmax::EnergyIota;

	// Do the final rounding.
	return efnmax::KCalToEnergy(logarithmic_part + strain_part + asymmetry_part);
}

efnmax::energy_t NNModelEFN2::MLAverageAsymmetryCost() const {
	return dt->mlasym;
}

efnmax::energy_t NNModelEFN2::MLStrainCost() const {
	return dt->strain;
}
