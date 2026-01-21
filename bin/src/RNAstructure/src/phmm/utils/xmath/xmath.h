#ifndef _XMATH_
#define _XMATH_

// Defines the double precision functions for linear and logarithmic math modes.
// This is an interface over xlog_math and linear_math functions. 
// The functions should make use of this class to do mathematical 
// operations in different modes at same time.
// The price to pay in this case is that the callback's will not be available for inlining
// which would most probably be substantial in case of linear computations since the type of
// computations are determined in run-time.
enum{MATH_MOD_LOG, MATH_MOD_LIN};

class t_xmath
{
public:
	t_xmath(int _math_mode);
	~t_xmath();

	int math_mode;

	// Callbacks for operations
	double (*sum)(double num1, double num2);
	double (*sub)(double num1, double num2);
	double (*div)(double num1, double num2);
	double (*mul)(double num1, double num2);
	double (*pow)(double base, double lin_exp);
	double (*max)(double num1, double num2);

	// Comparison operations
	bool (*compare)(double num1, double num2);
	bool (*geq)(double num1, double num2);
	bool (*gt)(double num1, double num2);
	
	// Compute the machine epsilon.
	double (*precision)();
};

#endif // _XMATH_