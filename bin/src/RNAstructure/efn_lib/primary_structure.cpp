//
// Created by Max Ward on 27 March 2017.
//

#include "primary_structure.hpp"

using namespace std;

efnmax::Base efnmax::CharToBase(char c) {
	// Note that we don't toupper(c) as a lower case
	// indicates a non-pairing nucleotide.
	switch (c) {
		case 'A':
			return A;
		case 'T':
			return U;
		case 'U':
			return U;
		case 'G':
			return G;
		case 'C':
			return C;
		default:
			break;
	}
	return NUMBASES; // Invalid base.
}

char efnmax::BaseToChar(efnmax::Base b) {
	switch (b) {
		case A:
			return 'A';
		case U:
			return 'U';
		case G:
			return 'G';
		case C:
			return 'C';
		case NUMBASES:
			return 'X'; // Special case for invalid base.
	}
	return 'X'; // this is not possible
}

efnmax::PrimeStructure efnmax::StringToPrimary(const string &seq) {
	efnmax::PrimeStructure bases;
	bases.reserve(seq.size());
	for (int i = 0; i < seq.size(); ++i) {
		char c = seq[i];
		bases.push_back(CharToBase(c));
	}
	return bases;
}

std::string efnmax::PrimaryToString(const PrimeStructure &primary) {
	string str;
	for (int i = 0; i < primary.size(); ++i) {
		Base b = primary[i];
		str.push_back(BaseToChar(b));
	}
	return str;
}
