#include <iostream>
#include <cmath>
#include <cfloat>
#include <limits>
#include <string>

const double LOG_OF_ZERO = -log(DBL_MAX)*1000;

double xlog(const double& d) {
  return d==0.0 ? LOG_OF_ZERO : log(d);
}
double xexp(const double& x) {
  return x<=LOG_OF_ZERO ? 0.0 : exp(x);
}

struct dbl {
	double val;
  dbl() { this->val = LOG_OF_ZERO; }
  dbl(const double& v) { this->val = xlog(v); }
	dbl(const double& v, const bool& b) { this->val = v; }
	double out() { return xexp(val); }
};

dbl xinit(const double& v) {
  dbl out;
  out.val = v;
  return out;
}

int main(int argc, char** argv) {
  if (argc < 2) { std::cerr << "Missing Arg1\n" << std::endl; return 1; }
  double v = std::stod(argv[1]);
  const char* name;
  #if D2
  name="D2";
  dbl a(xlog(v), true);
  #elif D3
  name="D3";
  dbl a;
  a = xinit(xlog(v));
  #else
  name="D1";
  dbl a = v;
  #endif

  std::cout << name << ":\t" << "a: " << a.val << "\t" << "A: " << a.out() << std::endl;
  return 0;
}
