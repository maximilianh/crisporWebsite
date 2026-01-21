//
// Created by Max Ward on 27 March 2017.
//

/**
 * Contains types and functions for secondary structure.
 * We define secondary structure as the result of base pairings,
 * but excluding pseudoknots for such pairings.
 */

#ifndef EFNMAX_SECONDARY_STRUCTURE_HPP
#define EFNMAX_SECONDARY_STRUCTURE_HPP

#include <vector>
#include <memory>
#include <string>
#include <cassert>

#include "primary_structure.hpp"

namespace efnmax {

/**
 * Matching represents a secondary structure.
 * The index i in the matching corresponds to the base i is bonded with.
 * In the case that i is unbonded, it will be paired with itself.
 * Another way of thinking about this that that a Matching is a permutation of integers [0, matching.size()).
 * This gives Matching some nice mathematical properties. It is closed under the "find paired" operator.
 * In addition, this sets up an involution: matching[matching[i]] == i.
 */
typedef std::vector<int> Matching;

/**
 * Creates and returns Matching of a particular size in which every base is unpaired.
 */
Matching EmptyMatching(unsigned size);

/**
 * Represents the paired nucleotides i and j. Assumed to be 5' to 3' so i<j.
 */
struct BondPair {
	// i and j are the nucleotides. Assumes i < j.
	int i, j;
	BondPair(int _i, int _j)
		: i(_i), j(_j) {
		assert(i < j);
	}
	/**
	 * A BondPair constructed this way represents "no bond"
	 */
	BondPair()
		: i(0), j(0) {}
	/**
	 * Returns true if i,j < other.j.
	 */
	bool operator<(const BondPair &other) const;
};
Matching DotBracketToMatching(const std::string &db);
std::string MatchingToDotBracket(const Matching &match);
/**
 * Whether two bases are a valid canonical bond pairing.
 */
bool ValidPair(Base a, Base b);
Matching BondPairsToMatching(const std::vector<BondPair> &bonds, unsigned sz);
std::vector<BondPair> MatchingToBondPairs(const Matching &match);
bool WatsonCrick(Base a, Base b);
bool ContainsPseudoknot(const Matching &match);

/**
 * Returns true iff the pair (i,j) has to be a lonley pair by nature of its neighbours.
 */
bool MustBeLonelyPair(const PrimeStructure &rna, int i, int j, int min_hairpin_unpaired);
}

#endif //EFNMAX_SECONDARY_STRUCTURE_HPP
