//Jacob Mainzer 5/10/09
//.cpp file for the extended_double class

//Including the precompile header fixes a link issue in mfc, see Microsoft Article ID: 148652
#ifdef _WINDOWS_GUI
#include "../Windows_interface_deprecated/stdafx.h"
#endif

#include "extended_double.h"
#include "rna_library.h"




#include <cmath>
#include <iostream>

using namespace std;

double extended_double::cap = 1.0e307;  //the cap size is predefined as this number.  1e307 is the largest possible cap size based on the implementation



//addition operator
extended_double operator+(const extended_double &ext1, const extended_double &ext2) {

	
	extended_double ADD;

	//3 cases: (2) if both terms have been scaled, (1) if one term has been scaled, (0) if neither term has been scaled
	switch (ext1.var_cap + ext2.var_cap)
	{
	case 2:
		//acts like normal addition
		ADD.var = ext1.var + ext2.var;
		ADD.var_cap = 1;
		//checks to see if sum dropped below cap limit (so that it can be scaled back down)
		if (ADD.var < 1 && ADD.var > -1) {
			ADD.var = ADD.var * extended_double::cap;  //scaled number scaled down by multiplying by cap
			ADD.var_cap = 0;
		}
		break;
	case 1:
		//first term is zero if the first number is not scaled, second term is zero if the second term is not scaled
		ADD.var = (ext1.var + ext2.var/extended_double::cap)*ext1.var_cap + (ext1.var/extended_double::cap + ext2.var)*ext2.var_cap;
		if (ADD.var < 1 && ADD.var > -1) {  //checks if sum is below the cap limit
			ADD.var = ADD.var * extended_double::cap;
			ADD.var_cap = 0;
		} else ADD.var_cap = 1;
		break;
	case 0:  //no scaled terms, so can act just like normal double addition, except for a cap limit check
		ADD.var = ext1.var + ext2.var;
		ADD.var_cap = 0;
		double temp = ADD.var/extended_double::cap; //a scaled version of ADD
		if (temp >= 1.0 || temp <= -1.0) {  //checks if sum is above cap limit
			ADD.var = temp;
			ADD.var_cap = 1;
		}
		break;
	}

	return ADD;
}


extended_double operator+(const extended_double & ext1, const double & ext2) {

	extended_double temp = ext2;

	return ext1 + temp;
}

extended_double operator+(const extended_double & ext1, const long double & ext2) {

	extended_double temp = ext2;

	return ext1 + temp;
}

//substraction is just addition with the second term turned negative
extended_double operator -(const extended_double & ext1, const extended_double & ext2) {

	return ext1 + (-ext2);
}

extended_double operator -(const extended_double & ext1, const double & ext2) {

	extended_double temp = ext2;

	return ext1 - temp;
}

//turns extended_double negative
extended_double operator -(extended_double ext1) {
	extended_double temp;
	temp.var = -ext1.var;
	temp.var_cap = ext1.var_cap;
	//temp.cap = cap;

	return temp;
}

//multiplication
extended_double operator *(const extended_double & ext1, const extended_double & ext2) {
	
	extended_double MULT;
	//cout << cap;
	
	//same idea as with addition, (2) is both are scaled, (1) one term is scaled, (0) no terms scaled
	switch (ext1.var_cap + ext2.var_cap)
	{
	case 2: //no need to check scaling since MULT cannot go below cap limit in this case
		MULT.var = ext1.var * ext2.var * extended_double::cap; //extra scaling term must be multiplied on to keep correct units
		MULT.var_cap = 1;
		break;
	case 1: //checks for drop below cap limit
		MULT.var = ext1.var * ext2.var;
		if (MULT.var < 1 && MULT.var > -1) {
			MULT.var = extended_double::cap * MULT.var;
			MULT.var_cap = 0;
		} else MULT.var_cap = 1;
		break;
	case 0: //uses scaled version of MULT to avoid overflows 
		double temp = (ext1.var / extended_double::cap) * ext2.var;
		if (temp >= 1 || temp <= -1) {  //scales
			MULT.var = temp;
			MULT.var_cap = 1;
		} else {  //doesnt scale
			MULT.var = ext1.var * ext2.var; 
			MULT.var_cap = 0;
		}
		break;
	}
	return MULT;
}

extended_double operator *(const extended_double & ext1, const double & ext2) {

	extended_double temp = ext2;

	return ext1 * temp;
}

//Division
extended_double operator /(const extended_double & ext1, const extended_double & ext2) {

	extended_double DIV;

	//same as before
	switch (ext1.var_cap + ext2.var_cap) {
		case 2:
			DIV.var = ext1.var / ext2.var;
			if (DIV.var < 1 && DIV.var > -1) {  //cap limit check
				DIV.var = DIV.var * extended_double::cap;
				DIV.var_cap = 0;
			} else DIV.var_cap = 1;
			break;
		case 1:  //first term is zero when first number is not scaled, second term is zero when second number is not scaled 
			DIV.var = ext1.var / ext2.var * ext1.var_cap + ((ext1.var / ext2.var) / extended_double::cap) * ext2.var_cap;
			if (DIV.var < 1 && DIV.var > -1) {
				if (ext2.var_cap == 0) DIV.var = DIV.var * extended_double::cap;
				DIV.var_cap = 0;
			} else DIV.var_cap = 1;
			break;	
		case 0:  //uses scaled version to catch overflows
			double temp = (ext1.var / extended_double::cap) / ext2.var;
			if (temp >=1 || temp <= -1) {
				DIV.var = temp;
				DIV.var_cap = 1;
			} else {
				DIV.var = ext1.var / ext2.var;
				DIV.var_cap = 0;
			}
			break;
	}
	return DIV;
}


extended_double operator /(const extended_double & ext1, const double & ext2) {

	extended_double temp = ext2;

	return ext1 / temp;
}

extended_double operator+(const double &ext1, const extended_double &ext2) {
	return ext2 + ext1;
}

extended_double operator-(const double &ext1, const extended_double &ext2) {
	extended_double temp = ext1;
	return temp  - ext2;
}

extended_double operator*(const double &ext1, const extended_double &ext2) {
	return ext2*ext1;
}

extended_double operator/(const double &ext1, const extended_double &ext2) {
	extended_double temp = ext1;
	return temp/ext2;
}





//internal write function
void extended_double::write() {
	if (var_cap == 0) cout << var;
	else {	
		int exp = (int)(floor(log10(fabs(var))));
		double mantissa = var*pow(10.0,-exp);
		exp += var_cap*(int)log10(cap);
		cout << mantissa << "e " << exp;
	}
}


bool operator<(const extended_double &ext1, const double &ext2) {
	if (ext1.var_cap == 0) return ext1.var < ext2;
	else return ext1.var < 0;  //if ext1 is scaled, then the only way it can be smaller than a double is if ext1 is negative.
				   //this logic is used repeatedly
}

bool operator<(const extended_double &ext1, const extended_double &ext2) {
	switch (ext1.var_cap + ext2.var_cap) {
		case 2: return ext1.var < ext2.var;
		case 1: if (ext1.var_cap == 1) return ext1.var < 0;
				else return ext2.var > 0;
		case 0: return ext1.var < ext2.var;
	}
	return false; //junk statement...never reached
}

bool operator<(const double &ext1, const extended_double &ext2) {
	if (ext2.var_cap == 0) return ext1 < ext2.var;
	else return ext2.var > 0;
}

bool operator==(const extended_double &ext1, const double &ext2) {
	if (ext1.var_cap == 0) return ext1.var == ext2;
	else return false;
}
bool operator==(const extended_double &ext1, const extended_double &ext2) {
	return (ext1.var_cap == ext2.var_cap)&&(ext1.var == ext2.var);
}
bool operator==(const double &ext1, const extended_double &ext2) {
	return  ext2 == ext1;
}


bool operator>(const extended_double &ext1, const double &ext2){
	return ext2 < ext1;
}

bool operator>(const extended_double &ext1, const extended_double &ext2){
	return ext2 < ext1;
}
bool operator>(const double &ext1, const extended_double &ext2){
	return ext2 < ext1;
}
bool operator<=(const extended_double &ext1, const double &ext2){
	return !(ext1 > ext2);
}
bool operator<=(const extended_double &ext1, const extended_double &ext2){
	return !(ext1 > ext2);
}
bool operator<=(const double &ext1, const extended_double &ext2){
	return !(ext1 > ext2);
}
bool operator>=(const extended_double &ext1, const double &ext2){
	return ext2 <= ext1;
}
bool operator>=(const extended_double &ext1, const extended_double &ext2){
	return ext2 <= ext1;
}
bool operator>=(const double &ext1, const extended_double &ext2){
	return ext2 <= ext1;
}

bool operator>(const extended_double &ext1, const int &ext2) {
	if (ext1.var_cap == 0) return ext1.var > ext2;
	else return ext1.var > 0;
}

bool operator>(const extended_double &ext1, const float &ext2) {
	if (ext1.var_cap == 0) return ext1.var > ext2;
	else return ext1.var > 0;
}

//lets you cout an extended double (ie cout << ext_double)
ostream& operator <<(ostream & out, const extended_double & ext){
	if (ext.var_cap == 0) out << ext.var;
	else {
		int exp = (int)(floor(log10(fabs(ext.var))));
		double mantissa = ext.var*pow(10.0,-exp);
		exp += ext.var_cap*(int)log10(extended_double::cap);
		out << mantissa << "e " << exp;
		
	}
	return out;
}
	



//overloaded read and write functions
void write(std::ofstream *out, extended_double *i) {
	//if (i->var_cap == 0) {
	//	out->write((char *) &(i->var), sizeof(i->var));
	//} else {
	//	int exp = (int)(floor(log10(i->var)) + i->var_cap*log10(i->cap));
	//	double mantissa = i->var*pow(10.0,-exp);
	//	out->write((char *) &(mantissa), sizeof(mantissa));
	//	out->write((char *) &(exp), sizeof(exp));
	//}
	write(out,&(i->var_cap));
	write(out,&(i->var));
}

//not sure this works
void read(std::ifstream *out, extended_double *i) {
	//if (i->var_cap == 0) {
    //            out->read((char *) &(i->var), sizeof(i->var));
    //    } else {
	//	int exp = (int)(floor(log10(i->var)) + i->var_cap*log10(i->cap));
	//	double mantissa = i->var*pow(10.0,-exp);
	//	out->read((char *) &(mantissa), sizeof(mantissa));
	//	out->read((char *) &(exp), sizeof(exp));
	//}
	read(out,&(i->var_cap));
	read(out,&(i->var));

}

//overloaded pow function
extended_double pow(const extended_double &base, const int &power) {
	
	if (base.var_cap == 0) return pow(base.var, power);

	extended_double temp = base;

	for (int i = 2; i <= power; i++){
		temp = temp * base;
	}
	return temp;
}

extended_double pow(const extended_double &base, const double &power) {
	
	if (base.var_cap == 0) return pow(base.var, power);

	extended_double temp = base;
	int tmp_power = (int) power;

	for (int i = 2; i <= tmp_power; i++){
		temp = temp * base;
	}
	return temp;
}

//overloaded log10 function
double log10(const extended_double &ext) {

	double temp = log10(ext.var) + ext.var_cap*log10(extended_double::cap);
	
	return temp;
}
//overloaded exp function
extended_double exp(const extended_double &ext) {

	extended_double temp = exp(ext.var);
	if (ext.var_cap == 1) temp = pow(temp, extended_double::cap); //always overflows
	return temp;
}
//overloaded log function
double log(const extended_double &ext) {
	double temp = log(ext.var) + ext.var_cap*log(extended_double::cap);
	
	return temp;
}
