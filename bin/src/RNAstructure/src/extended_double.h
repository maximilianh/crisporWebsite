

//Jacob Mainzer 5/10/09
//Extended double class.  Uses a scaled double to extend the maximum limit to ~1e615
//uses a scaling method and a cap so that when the double becomes above the cap size, the double is divided
//by the size of the cap and marker variable is changed to denote a scaling.

//Add check so that the header is included only once:
#ifndef ED_H
#define ED_H


#include <fstream>
#include <iostream>
#include <cmath>


//using namespace std;

class extended_double {

public:

	double var;  //double 	
	short var_cap;  //is 0 if not scaled, 1 if scaled by the cap size	
	static double cap; //the cap size
	

	//constructors
	extended_double(){
	}
	extended_double(const extended_double& ext) {
		var = ext.var;
		var_cap = ext.var_cap;
	}

	extended_double(const double& ext) {
		var = ext;
		var_cap = 0;
	}

	//operators using various combinations of types of variables
	extended_double& operator=(const extended_double &ext2);
	extended_double& operator=(const double &ext2);
	extended_double& operator=(const long double &ext2);
	friend extended_double operator+(const extended_double &ext1, const extended_double &ext2);
	friend extended_double operator-(const extended_double &ext1, const extended_double &ext2);
	friend extended_double operator*(const extended_double &ext1, const extended_double &ext2);
	friend extended_double operator/(const extended_double &ext1, const extended_double &ext2);
	friend extended_double operator+(const extended_double &ext1, const double &ext2);
	friend extended_double operator+(const extended_double &ext1, const long double &ext2);
	friend extended_double operator-(const extended_double &ext1, const double &ext2);
	friend extended_double operator*(const extended_double &ext1, const double &ext2);
	friend extended_double operator/(const extended_double &ext1, const double &ext2);
	friend extended_double operator+(const double &ext1, const extended_double &ext2);
	friend extended_double operator-(const double &ext1, const extended_double &ext2);
	friend extended_double operator*(const double &ext1, const extended_double &ext2);
	friend extended_double operator/(const double &ext1, const extended_double &ext2);
	friend extended_double operator-(extended_double ext1);
	const extended_double& operator+=(const extended_double &ext2);
	const extended_double& operator+=(const double &ext2);
	friend bool operator<(const extended_double &ext1, const double &ext2);
	friend bool operator<(const extended_double &ext1, const extended_double &ext2);
	friend bool operator<(const double &ext1, const extended_double &ext2);

	friend bool operator>(const extended_double &ext1, const double &ext2);
	friend bool operator>(const extended_double &ext1, const extended_double &ext2);
	friend bool operator>(const double &ext1, const extended_double &ext2);
	friend bool operator>(const extended_double &ext1, const int &ext2);
	friend bool operator>(const extended_double &ext1, const float &ext2);
	
	friend bool operator<=(const extended_double &ext1, const double &ext2);
	friend bool operator<=(const extended_double &ext1, const extended_double &ext2);
	friend bool operator<=(const double &ext1, const extended_double &ext2);
	
	friend bool operator>=(const extended_double &ext1, const double &ext2);
	friend bool operator>=(const extended_double &ext1, const extended_double &ext2);
	friend bool operator>=(const double &ext1, const extended_double &ext2);
	
	friend bool operator==(const extended_double &ext1, const double &ext2);
	friend bool operator==(const extended_double &ext1, const extended_double &ext2);
	friend bool operator==(const double &ext1, const extended_double &ext2);
	friend std::ostream& operator <<(std::ostream & out, const extended_double & ext);

	void write();
	
	//type cast ##### will try to return the number of the extended double in a double form, so it could overflow if the number is too large####
	operator double(){
		//cout << "double conversion: " << var << " " << var_cap << endl;
		return var_cap*var*extended_double::cap + (1 - var_cap)*var;
	}

private:
	

};

//overloaded functions
void write(std::ofstream *out, extended_double *i);
void read(std::ifstream *out, extended_double *i);
extended_double pow(const extended_double &base, const int &power);
extended_double pow(const extended_double &base, const double &power);
double log10(const extended_double &ext);
double log(const extended_double &ext);
extended_double exp(const extended_double &ext);



inline extended_double& extended_double::operator=(const extended_double &ext2) {

	var = ext2.var;
	var_cap = ext2.var_cap;

	//cap = 1e307;

	return *this;
}

inline extended_double& extended_double::operator=(const double &ext2) {

	var = ext2;
	var_cap = 0;

	

	return *this;
}

inline extended_double& extended_double::operator=(const long double &ext2) {

	var = ext2;
	var_cap = 0;

	

	return *this;
}




inline const extended_double& extended_double::operator +=(const extended_double &ext2) {

	*this = *this + ext2;

	return *this;
}

inline const extended_double& extended_double::operator +=(const double &ext2) {

	extended_double temp;
	temp = ext2;

	*this = *this + temp;

	return *this;

}



#endif //ED_H


