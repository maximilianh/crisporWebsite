//
// Created by Max Ward on 27 March 2017.
//

/**
 * Contains functions and types related to RNA primary structure.
 * The primary structure of an RNA is the nucleotide sequence.
 */

#ifndef EFNMAX_PRIMARY_STRUCTURE_HPP
#define EFNMAX_PRIMARY_STRUCTURE_HPP

#include <vector>
#include <string>

namespace efnmax {
/**
  * Represents a base on a nucleotide.
  */
enum Base { A = 0, U, G, C, NUMBASES };
/**
  * The primary structure of an RNA. The sequence of bases itself.
  */
typedef std::vector<Base> PrimeStructure;
/**
  * Converts a given character into the corresponding Base.
  */
Base CharToBase(char c);
/**
  * Converts a base into the corresponding character.
  */
char BaseToChar(Base b);
/**
  * Parses a string into a primary structure.
  */
PrimeStructure StringToPrimary(const std::string &seq);
/**
  * Represents a primary structure as a string..
  */
std::string PrimaryToString(const PrimeStructure &primary);
}

#endif //EFNMAX_PRIMARY_STRUCTURE_HPP
