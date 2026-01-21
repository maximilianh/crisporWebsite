#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <limits>

const double LOG_OF_ZERO = -std::numeric_limits<double>::infinity();

double lsum(const double log1, const double log2) {
	return log1 + log(1 + exp(log2-log1));
}
double lsum1p(const double log1, const double log2) {
	return log1 + log1p(exp(log2-log1));
}
double lsum1p_min(const double log1, const double log2) {
	return log1>log2?log1+log1p(exp(log2-log1)):log2+log1p(exp(log1-log2));
}

double xlog_mul(double log1, double log2) {
	return log1+log2;
}

// Returns 0 if log1 is 0 no matter what log2 is.
double xlog_div(double log1, double log2) {
	return log1-log2;
}

inline double xlog_sum(double a, double b) {
	return a>b? a+log1p(exp(b-a)):  b+log1p(exp(a-b)); // log1pexp performs a check to see if its argument is <= LOG_OF_ZERO
}

// Computes log(exp(a)-exp(b))
inline double xlog_sub(double a, double b) { // (previously xlog_sub)
	return a+log1p(-exp(b-a));  // note that (b-a) is always < 0, which is required because we are computing log(1-exp(b-a))
}

void compare(double(*f)(const double, const double), const double a, const double b, const double sum) {
	double ra=f(a,b),      rb=f(b,a);
	double sa=exp(ra),     sb=exp(rb);
	double da=abs(sum-sa), db=abs(sum-sb);
	std::cout << "(a,b)="<<ra<<"\tsum="<<sa<<"\tdelta="<<da<<"\n";
	std::cout << "(b,a)="<<rb<<"\tsum="<<sb<<"\tdelta="<<db<<"\n";
	std::cout << (da==db?"Same":da<db?"(a,b) better":"(b,a) better")<<" (diff: "<<abs(da-db)<<")\n";
}

void eval(const char* const name, double(*f)(const double, const double), const double A, const double B) {
	double a=log(A), b=log(B), c;
	std::cout<<"Test "<<name<<"\n";
	std::cout << "a="<<a<<"\tb="<<b<<"\n";
	c = f(a,b);
	std::cout << "c="<<c<<"\tC="<<exp(c)<<"\n";
	std::cout.flush();
}

int main(int argc, char** argv) {
	if (argc < 3) { std::cerr<<"Requires two args!"<<std::endl; return 1; }
	double A = atof(argv[1]), B = atof(argv[2]);
	double a=log(A), b=log(B);
	int runs = 100;
	if (argc >=4) runs = atoi(argv[3]);

	std::cout << std::setprecision(18);
	std::cout << "A="<<A<<"\tlog(A)="<<a<<"\n";
	std::cout << "B="<<B<<"\tlog(B)="<<b<<"\n";
	double sum = A+B;
	std::cout << "SUM="<< sum<<"\n";
	
	std::cout << "----lsum----\n";
	compare(lsum,a,b,sum);

	std::cout << "----lsum1p----\n";
	compare(lsum1p,a,b,sum);

	std::cout << "----lsum1p_min----\n";
	compare(lsum1p_min,a,b,sum);	


	eval("0+1",xlog_sum,0,1);
	eval("0*1",xlog_mul,0,1);

	eval("1+0",xlog_sum,1,0);
	eval("1*0",xlog_mul,1,0);

	eval("0/1",xlog_div,0,1);
	eval("1-0",xlog_sub,1,0);


	eval("0+1E-10",xlog_sum,0,1E-10);
	eval("0*1E-10",xlog_mul,0,1E-10);

	eval("1E-10+0",xlog_sum,1E-10,0);
	eval("1E-10*0",xlog_mul,1E-10,0);

	eval("0/1E-10",xlog_div,0,1E-10);
	eval("1E-10-0",xlog_sub,1E-10,0);
	
	eval("1-1",xlog_sub,1,1);
	eval("1000-1000",xlog_sub,1000,1000);
	eval("1E-10-1E-10",xlog_sub,1E-10,1E-10);

	// double a1=a, a2=a, A1=A, A2=A;
	// double b1=b, b2=b, B1=B, B2=B;
	// double c=b, C=B;
	// for (int i = 0; i < runs; i++) {
	// 	a1 = lsum(a1,b);
	// 	a2 = lsum(b, a2);
	// 	A1 = A1 + B;
	// 	A2 = B + A2;

	// 	b2 = lsum(a, b2);
	// 	b1 = lsum(b1,a);
	// 	B2 = A + B2;
	// 	B1 = B1 + A;
	// }
	// for (int i = 0; i < (runs-1); i++) {
	// 	c = lsum(c, b);
	// 	C = C+B;
	// }
	// double ABB=(A+(runs*B)), BAA=(B+(runs*A));

	// std::cout << "\nRepeated Sums\n";
	// std::cout << "A+("<<runs<<"*B)=" << ABB << "\n";
	// std::cout << "B+("<<runs<<"*A)=" << BAA << "\n";

	// std::cout << "\nStandard Additions\n";
	// std::cout << "A1=" << A1 << "\tdelta=" << abs(ABB-A1) << "\n";
	// std::cout << "A2=" << A2 << "\tdelta=" << abs(ABB-A2) << "\n";
	// std::cout << "B1=" << B1 << "\tdelta=" << abs(BAA-B1) << "\n";
	// std::cout << "B2=" << B2 << "\tdelta=" << abs(BAA-B2) << "\n";

	// std::cout << "\nLog-base Additions\n";
	// std::cout << "a1=" << a1 << "\tVal=" << exp(a1) << "\tdelta=" << abs(ABB-exp(a1)) << "\n";
	// std::cout << "a2=" << a2 << "\tVal=" << exp(a2) << "\tdelta=" << abs(ABB-exp(a2)) << "\n";
	// std::cout << "b1=" << b1 << "\tVal=" << exp(b1) << "\tdelta=" << abs(BAA-exp(b1)) << "\n";
	// std::cout << "b2=" << b2 << "\tVal=" << exp(b2) << "\tdelta=" << abs(BAA-exp(b2)) << "\n";

	// std::cout << "\nAggregated Additions\n";
	// std::cout << "A+C=" << (A+C) << "\n";
	// std::cout << "a+c=" << lsum(a, c) << "\tVal=" << exp(lsum(a, c)) << "\tdelta=" << abs((A+C)-exp(lsum(a, c))) << "\n";

	std::cout.flush();

	return 0;
}