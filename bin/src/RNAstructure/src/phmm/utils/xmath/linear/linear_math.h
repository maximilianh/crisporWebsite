#ifndef _LINEAR_MATH_
#define _LINEAR_MATH_

//#define EPS (pow(0.1, 10))

//#define MAX(x, y) ((x)>(y)?(x):(y)) // MIN and MAX are now defined in phmm.h
//#define MIN(x, y) ((x)>(y)?(y):(x))

double lin_sum(double num1, double num2);
double lin_sub(double num1, double num2);
double lin_mul(double num1, double num2);
double lin_div(double num1, double num2);
double lin_pow(double num1, double num2);
double lin_max(double num1, double num2);

bool lin_compare(double num1, double num2);
bool lin_geq(double num1, double num2);
bool lin_gt(double val1, double val2);

double get_linear_comp_prec();

#endif // _LINEAR_MATH_


