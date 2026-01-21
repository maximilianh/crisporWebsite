#include "mainloop.h"
#include "sequence.h"
#include "NCM_parameters.h"
#include <vector>

std::vector<table_t> turbo_calculation(const std::vector<sequence>& seqs,
                                       const parameters<real_t>& p,
                                       const int iterations,
                                       const real_t gamma=0.6,
                                       const std::vector<constraints>& constr=std::vector<constraints>());
