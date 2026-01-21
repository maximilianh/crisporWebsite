//
// Created by Max Ward on 27 March 2017.
//

#include "nn_model.hpp"
#include "secondary_structure.hpp"

#include <cassert>
#include <limits>

using namespace std;

efnmax::energy_t NNModel::OneLoop(int i, int j) const {
	assert(i < j);
	assert(efnmax::ValidPair(rna[i], rna[j]));
	if (j - i - 1 < MIN_HAIRPIN_UNPAIRED)
		return MaxMFE();
	return erg3(i + 1, j + 1, struc, dt, 0);
}

efnmax::energy_t NNModel::TwoLoop(int i, int k, int l, int j) const {
	assert(i < j && k < l && i < k && l < j);
	if (k == i + 1 && l == j - 1)
		return erg1(i + 1, j + 1, k + 1, l + 1, struc, dt);
	else
		return erg2(i + 1, j + 1, k + 1, l + 1, struc, dt, 0, 0);
}

efnmax::energy_t NNModel::Branch(int i, int j) const {
	return penalty(i + 1, j + 1, struc, dt);
}


efnmax::energy_t NNModel::FlushCoax(int i, int j, int k, int l) const {
	assert(i < j && k < l && (k == j + 1 || k == i + 1 || l == j - 1));
	if (i < k && l < j) // i,j close a multi-loop.
		return erg1(i + 1, j + 1, k + 1, l + 1, struc, dt);
	return erg1(j + 1, i + 1, k + 1, l + 1, struc, dt);
}

efnmax::energy_t NNModel::MismatchCoax(int i, int j, int k, int l) const {
	assert(i < j && k < l && (abs(i - k) == 2 || abs(i - l) == 2 || abs(j - k) == 2 || abs(j - l) == 2));
	// i,j is mismatched stacked with k,l
	// Mismatch is off i,j
	i += 1;
	j += 1;
	k += 1;
	l += 1;
	if (i < k && l < j) { // i,j closes a multiloop
		if (k == i + 2) // (.(_)_.)
			return ergcoaxinterbases1(j, i, k, l, struc, dt);
		else // (._(_).)
			return ergcoaxinterbases2(k, l, j, i, struc, dt);
	} else if (k < i && j < l) { // k,l closes a multiloop
		if (i == k + 2) // (.(_)._)
			return ergcoaxinterbases2(l, k, i, j, struc, dt);
		else // (_.(_).)
			return ergcoaxinterbases1(i, j, l, k, struc, dt);
	} else if (j < k) { // .(_).(_)
		return ergcoaxinterbases1(i, j, k, l, struc, dt);
	} else { // (_).(_).
		return ergcoaxinterbases2(k, l, i, j, struc, dt);
	}

}

efnmax::energy_t NNModel::FiveDangle(int i, int j) const {
	assert(i < j && i > 0);
	return erg4(j + 1, i + 1, i, 2, struc, dt, false);
}

efnmax::energy_t NNModel::ClosingFiveDangle(int i, int j) const {
	assert(i < j);
	return erg4(i + 1, j + 1, j, 2, struc, dt, false);
}

efnmax::energy_t NNModel::ThreeDangle(int i, int j) const {
	assert(i < j && j + 1 < static_cast<int>(rna.size()));
	return erg4(j + 1, i + 1, j + 2, 1, struc, dt, false);
}

efnmax::energy_t NNModel::ClosingThreeDangle(int i, int j) const {
	assert(i < j);
	return erg4(i + 1, j + 1, i + 2, 1, struc, dt, false);
}

efnmax::energy_t NNModel::Mismatch(int i, int j) const {
	assert(i < j);
	return dt->tstkm[struc->numseq[j + 1]][struc->numseq[i + 1]][struc->numseq[j + 2]][struc->numseq[i]];
}

efnmax::energy_t NNModel::ClosingMismatch(int i, int j) const {
	assert(i < j);
	return dt->tstkm[struc->numseq[i + 1]][struc->numseq[j + 1]][struc->numseq[i + 2]][struc->numseq[j]];
}

void NNModel::SetRNA(const std::string &primary) {
	this->rna = efnmax::StringToPrimary(primary);
	delete this->struc;
	this->struc = new structure();

	this->struc->allocate(static_cast<int>(primary.size()));
	for (int i = 1; i <= this->struc->GetSequenceLength(); ++i) {
		efnmax::Base nuc = efnmax::CharToBase(primary[i - 1]);
		if (nuc == efnmax::A) this->struc->numseq[i] = 1;
		else if (nuc == efnmax::C) this->struc->numseq[i] = 2;
		else if (nuc == efnmax::G) this->struc->numseq[i] = 3;
		else if (nuc == efnmax::U) this->struc->numseq[i] = 4;
		else this->struc->numseq[i] = 0;
		this->struc->nucs[i] = efnmax::BaseToChar(nuc);
		this->struc->hnumber[i] = static_cast<short>(i);
	}
}

efnmax::energy_t NNModel::MaxMFE() const {
	return numeric_limits<efnmax::energy_t>::max() / 3;
}