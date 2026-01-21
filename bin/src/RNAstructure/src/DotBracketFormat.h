#ifndef DOT_BRACKET_FORMAT_H
#define DOT_BRACKET_FORMAT_H

#include <string>
#include "common_utils.h"

// ************************************************
//  NOTE: These definitions are used by the writedotbracket function in the structure class.
//  These could be included in structure.h, have been moved to this separate header 
//  to avoid code bloat in structure.h
// ************************************************

// Bit definitions for dot-bracket notation output format
#define DBN_BIT_SEQ_LABEL 0x1    // write the sequence label instead of the sub-structure label in the title-line
#define DBN_BIT_SIDE_TITLES 0x2  // write sub-structure titles on the side of the dot-bracket structure line
#define DBN_BIT_MULTI_TITLE 0x4  // write a ">title" line for each sub-structure
#define DBN_BIT_MULTI_SEQ 0x8    // write a sequence line for each sub-structure

/**
  *  DotBracketFormat -- an enum that describes how multiple structures should be written. 
  *    A dot-bracket file begins with a title line, followed by a sequence line, and then 
  *    a structure line (see definitions below). If a dot-bracket file contains MULTIPLE structures, 
  *    there is a choice about how to write subsequent sub-structures.
  *
  *  Style 1 (DBN_FMT_SINGLE_TITLE)
  *    Subsequent structures are listed only as additional structure lines (with no titles).
  *    This format results in a loss of information if each sub-structure has a unique label/title.
  *
  *  Style 2 (DBN_FMT_SIDE_TITLES)
  *    Subsequent structures are listed as additional structure lines with a side comment -- i.e. the 
  *    structure line is followed by a tab character and then a comment/title (all on the same line).
  *    This format preserves sub-structure labels (if present).
  *
  *  Style 3 (DBN_FMT_MULTI_TITLE)
  *    Subsequent structures are written starting with an individual title line followed by a structure line. 
  *    I.e. -- just like the first structure except that the sequence line is NOT repeated.
  *
  *  Style 2 (DBN_FMT_MULTI_TITLE_AND_SEQ)
  *    Subsequent structures are listed the same as the first structure -- i.e. a title line followed by a 
  *    sequence line and then a structure line.
  *
  *  #### Definitions -- components of a dot-bracket file ####
  *    Title Line     -- contains '>' followed by title of structure:   >ENERGY = -33.6  RA7680
  *    Sequence Line  -- contains nucleobase sequence. e.g. :           GGCCUAUGGCC
  *    Structure Line -- basepair information in dot-bracket notation:  ((((...))))
  *    Structure Line with side comment/title:                          ((((...))))  ENERGY = -33.6
  *
  */
enum DotBracketFormat {
	//! Subsequent structures are written with structure-lines only.
	DBN_FMT_SINGLE_TITLE = DBN_BIT_SEQ_LABEL,      
	//! Stucture lines have side-titles
	DBN_FMT_SIDE_TITLES  = DBN_BIT_SEQ_LABEL|DBN_BIT_SIDE_TITLES,  
	//! Subsequent structures have >title lines.
	DBN_FMT_MULTI_TITLE  = DBN_BIT_MULTI_TITLE,      
	//! Subsequent structures have >title lines and repeated sequences.
	DBN_FMT_MULTI_TITLE_AND_SEQ  = DBN_BIT_MULTI_TITLE|DBN_BIT_MULTI_SEQ 
};

inline DotBracketFormat parseDotBracketFormat(const std::string &format) {
    const std::string fmt = toLower(format);
    if (fmt=="simple"||fmt=="1") return DBN_FMT_SINGLE_TITLE;
    if (fmt=="side"  ||fmt=="2") return DBN_FMT_SIDE_TITLES;
    if (fmt=="multi" ||fmt=="3") return DBN_FMT_MULTI_TITLE;
    if (fmt=="full"  ||fmt=="4") return DBN_FMT_MULTI_TITLE_AND_SEQ;
    return static_cast<DotBracketFormat>(0); // indicates error
}

#endif //DOT_BRACKET_FORMAT_H