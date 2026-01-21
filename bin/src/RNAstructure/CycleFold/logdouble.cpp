#include "logdouble.h"
#include <cmath>
#include <cfloat>
#include <cassert>

const double LOG_ALMOST_ZERO = -log(DBL_MAX);
const double LOG_ALMOST_INF = log(DBL_MAX);

inline double log1pexp(double x)
{
	return x < LOG_ALMOST_ZERO? 0.0 : log1p(exp(x));
}

inline double log_sum(double a, double b)
{
    if(a<=LOG_ALMOST_ZERO) return b;
    if(b<=LOG_ALMOST_ZERO) return a;
	return a>b? a+log1pexp(b-a):  b+log1pexp(a-b);
}

inline double log_diff(double a, double b)
{
	assert(a>=b);
	return a + log1p(-exp(b-a));
}

log_double::log_double() : val(-DBL_MAX)
{}

log_double::log_double(const double& a) : val(log(a)){}

log_double::log_double(const int& a) : val(log((double) a)){}

log_double::log_double(const log_double& a) : val(a.val){}

log_double::operator double(){
	return exp(val);
}

log_double::operator float(){
	return exp(val);
}

double log_double::asDouble() const{
	if (val>LOG_ALMOST_INF){
		std::cout<<"possible overflow converting from log\n";
	}
	return exp(val);
}

const log_double& log_double::operator=(const log_double& rhs)
{
	val = rhs.val;
	return *this;
}

log_double log_double::operator+(const log_double& rhs) const
{
	log_double ret(*this);
	ret += rhs;
	return ret;
}

log_double log_double::operator-(const log_double& rhs) const
{
	log_double ret(*this);
	ret -= rhs;
	return ret;
}

double log_double::operator-() const{
    return -val;
}

log_double log_double::operator*(const log_double& rhs) const
{
	log_double ret(*this);
	ret *= rhs;
	return ret;
}

log_double log_double::operator/(const log_double& rhs) const
{
	log_double ret(*this);
	ret /= rhs;
	return ret;
}

log_double log_double::operator^(const int& rhs) const{
    log_double ret = 1.0;
    for(int i=0;i<rhs;i++){
        ret *= *this;
    }
    return ret;
}

const log_double& log_double::operator+=(const log_double& rhs)
{
	this->val = log_sum(this->val,rhs.val);
	return *this;
}

const log_double& log_double::operator-=(const log_double& rhs)
{
	this->val = log_diff(this->val,rhs.val);
	return *this;
}

const log_double& log_double::operator*=(const log_double& rhs)
{
	this->val = this->val + rhs.val;
	return *this;
}

const log_double& log_double::operator/=(const log_double& rhs)
{
	this->val = this->val - rhs.val;
	return *this;
}

bool log_double::operator==(const log_double& rhs) const
{
	return this->val == rhs.val;
}

bool log_double::operator>(const log_double& rhs) const
{
	return this->val > rhs.val;
}

bool log_double::operator<(const log_double& rhs) const
{
	return this->val < rhs.val;
}

bool log_double::operator>=(const log_double& rhs) const
{
	return this->val >= rhs.val;
}

bool log_double::operator<=(const log_double& rhs) const
{
	return this->val <= rhs.val;
}

std::ostream& operator<<(std::ostream& os, const log_double& ld)
{
  return os << ld.asDouble();
}

//arithmetic with doubles: always promote the double to a log_double
//then return the result using log_double operation
log_double log_double::operator+(const double& rhs) const
{
    return *this + log_double(rhs);
}

log_double log_double::operator-(const double& rhs) const
{
    return *this - log_double(rhs);
}

log_double log_double::operator*(const double& rhs) const
{
    return *this * log_double(rhs);
}

log_double log_double::operator/(const double& rhs) const
{
    return *this / log_double(rhs);
}

bool log_double::operator==(const double& rhs) const
{
    return *this == log_double(rhs);
}

bool log_double::operator>(const double& rhs) const
{
    return *this > log_double(rhs);
}

bool log_double::operator<(const double& rhs) const
{
    return *this < log_double(rhs);
}

bool log_double::operator>=(const double& rhs) const
{
    return *this >= log_double(rhs);
}

bool log_double::operator<=(const double& rhs) const
{
    return *this <= log_double(rhs);
}

log_double operator+(const double& lhs, const log_double& rhs)
{
    return log_double(lhs) + rhs;
}

log_double operator-(const double& lhs, const log_double& rhs)
{
    return log_double(lhs) - rhs;
}

log_double operator*(const double& lhs, const log_double& rhs)
{
    return log_double(lhs) * rhs;
}

log_double operator/(const double& lhs, const log_double& rhs)
{
    return log_double(lhs) / rhs;
}

bool operator==(const double& lhs, const log_double& rhs)
{
    return log_double(lhs) == rhs;
}

bool operator>(const double& lhs, const log_double& rhs)
{
    return log_double(lhs) > rhs;
}

bool operator<(const double& lhs, const log_double& rhs)
{
    return log_double(lhs) < rhs;
}

bool operator>=(const double& lhs, const log_double& rhs)
{
    return log_double(lhs) >= rhs;
}

bool operator<=(const double& lhs, const log_double& rhs)
{
    return log_double(lhs) <= rhs;
}

bool log_double::operator==(const int& rhs) const
{
    return *this == log_double(rhs);
}

bool log_double::operator>(const int& rhs) const
{
    return *this > log_double(rhs);
}

bool log_double::operator<(const int& rhs) const
{
    return *this < log_double(rhs);
}

bool log_double::operator>=(const int& rhs) const
{
    return *this >= log_double(rhs);
}

bool log_double::operator<=(const int& rhs) const
{
    return *this <= log_double(rhs);
}

bool operator==(const int& lhs, const log_double& rhs)
{
    return log_double(lhs) == rhs;
}

bool operator>(const int& lhs, const log_double& rhs)
{
    return log_double(lhs) > rhs;
}

bool operator<(const int& lhs, const log_double& rhs)
{
    return log_double(lhs) < rhs;
}

bool operator>=(const int& lhs, const log_double& rhs)
{
    return log_double(lhs) >= rhs;
}

bool operator<=(const int& lhs, const log_double& rhs)
{
    return log_double(lhs) <= rhs;
}


log_double log(const log_double& ld)
{
    return log_double(log(ld.val));
}

log_double exp(const log_double& ld)
{
    return log_double(exp(ld.val));
}

log_double pow(const log_double& ld, const int& exponent)
{
    double result = pow(ld.asDouble(), (double) exponent);
    return log_double(result);
}

log_double pow(const log_double& ld, const log_double& exponent)
{
    double result = pow(ld.asDouble(), exponent.asDouble());
    return log_double(result);
}

void read(std::ifstream *out, log_double* i)
{
    out->read((char *) i, sizeof(*i));
}
void write(std::ofstream *out, log_double* i)
{
    out->write((char *) i, sizeof(*i));
}
