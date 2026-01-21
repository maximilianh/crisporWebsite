#ifndef __LOGDOUBLE__
#define __LOGDOUBLE__
#include <iostream>
#include <fstream>
//this class allows convenient manipulation of real values in log space
//the operators follow the laws of logarithms, so for example multiplication
//will add the underlying values
class log_double
{
public:
	log_double();
	log_double(const double& a);
	log_double(const int& a);
	log_double(const log_double& a);
	double asDouble() const; //get the exponential of val in a double
	const log_double& operator=(const log_double& rhs);
	log_double operator+(const log_double& rhs) const;
	log_double operator-(const log_double& rhs) const;
	log_double operator*(const log_double& rhs) const;
	log_double operator/(const log_double& rhs) const;
	double operator-() const;
	log_double operator^(const int& rhs) const;
	const log_double& operator+=(const log_double& rhs);
	const log_double& operator-=(const log_double& rhs);
	const log_double& operator*=(const log_double& rhs);
	const log_double& operator/=(const log_double& rhs);
	bool operator==(const log_double& rhs) const;
	bool operator>(const log_double& rhs) const;
	bool operator<(const log_double& rhs) const;
    bool operator>=(const log_double& rhs) const;
    bool operator<=(const log_double& rhs) const;
//#ifdef NON_TYPESAFE_LOGDOUBLE
    //use for compatibility with old RNAstructure code
    //would really be better to fix things so it's not necessay
//arithmetic operators
	log_double operator+(const double& rhs) const;
	log_double operator-(const double& rhs) const;
	log_double operator*(const double& rhs) const;
	log_double operator/(const double& rhs) const;
	friend log_double operator+(const double& lhs,const log_double& rhs);
	friend log_double operator-(const double& lhs,const log_double& rhs);
	friend log_double operator*(const double& lhs,const log_double& rhs);
	friend log_double operator/(const double& lhs,const log_double& rhs);
	const log_double& operator+=(const double& rhs);
	const log_double& operator-=(const double& rhs);
	const log_double& operator*=(const double& rhs);
	const log_double& operator/=(const double& rhs);
    //if we need these then something is terribly wrong
/*	friend const log_double& operator+=(const double& lhs, const log_double& rhs);
	friend const log_double& operator-=(const double& lhs, const log_double& rhs);
	friend const log_double& operator*=(const double& lhs, const log_double& rhs);
	friend const log_double& operator/=(const double& lhs, const log_double& rhs);
*/
//comparison operators
//compare to double
    bool operator==(const double& rhs) const;
	bool operator>(const double& rhs) const;
	bool operator<(const double& rhs) const;
    bool operator>=(const double& rhs) const;
    bool operator<=(const double& rhs) const;
	friend bool operator==(const double& lhs, const log_double& rhs);
	friend bool operator>(const double& lhs, const log_double& rhs);
	friend bool operator<(const double& lhs, const log_double& rhs);
    friend bool operator>=(const double& lhs, const log_double& rhs);
    friend bool operator<=(const double& lhs, const log_double& rhs);

//compare to int
    bool operator==(const int& rhs) const;
	bool operator>(const int& rhs) const;
	bool operator<(const int& rhs) const;
    bool operator>=(const int& rhs) const;
    bool operator<=(const int& rhs) const;
	friend bool operator==(const int& lhs, const log_double& rhs);
	friend bool operator>(const int& lhs, const log_double& rhs);
	friend bool operator<(const int& lhs, const log_double& rhs);
    friend bool operator>=(const int& lhs, const log_double& rhs);
    friend bool operator<=(const int& lhs, const log_double& rhs);

//#endif
	operator double();
	operator float();
	double val; //holds the log value

    friend log_double log(const log_double& ld);
    friend log_double exp(const log_double& ld);
    friend log_double pow(const log_double& ld, const int& exponent);
    friend log_double pow(const log_double& ld, const log_double& exponent);
    friend void read(std::ifstream *out, log_double* i);
    friend void write(std::ofstream *out, log_double* i);
};


typedef log_double real_t;

std::ostream& operator<<(std::ostream& os, const log_double& ld);

#endif
