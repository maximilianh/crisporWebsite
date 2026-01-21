//
// Created by Max Ward on 29 March 2017.
//

#include "nn_model_linear.hpp"

#include <cassert>


efnmax::energy_t NNModelLinear::MLClosure(unsigned branches, unsigned unpaired) const {
	assert(branches >= 3);
	return MLInitCost() + branches*MLBranchCost() + unpaired*MLUnpairedCost();
}

efnmax::energy_t NNModelLinear::MLInitCost() const {
	return dt->efn2a;
}

efnmax::energy_t NNModelLinear::MLBranchCost() const {
	return dt->efn2c;
}

efnmax::energy_t NNModelLinear::MLUnpairedCost() const {
	return dt->efn2b;
}