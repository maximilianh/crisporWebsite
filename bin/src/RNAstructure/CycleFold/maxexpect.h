#ifndef GENERATE_STRUCTURES_H
#define GENERATE_STRUCTURES_H
#include "mainloop.h"

table_t maxexpect_fill(const table_t& p, const std::vector<double>& q, double gamma);
std::vector<std::pair<int,int>> maxexpect_traceback(const table_t& M, const table_t& p, const std::vector<double>& q, const double gamma);
std::vector<std::pair<int,int>> maxexpect(const table_t& p, double gamma=1.0);
#endif
