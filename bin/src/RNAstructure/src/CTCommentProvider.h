#ifndef CT_COMMENT_PROVIDER_H
#define CT_COMMENT_PROVIDER_H

#include "defines.h"
#include <string>
#include <cstdio>

// ************************************************
//  NOTE: the implemenation for this class is included in structure.cpp
//  The only reason this separate header exists is to avoid code bloat in structure.h
// ************************************************

class structure; //forward declaration of structure class 

#define NO_CT_COMMENTS     CTComments::None
#define CT_ENERGY_COMMENTS CTComments::Energy

//!     This class allows a program to customize comments written for each structure in a CT or dot-bracket file.
//!     Many programs write "ENERGY = ..." (when GetEnergy() is not zero).
//!     However, a program like MaxExpect might write "SCORE = ..." (even if GetEnergy() IS zero).
//!     Similarly, Design might want to add comments like "NED = 0.02" which would require a different source of information than
//!     GetEnergy() -- or at the very least, different formatting/precision. 
//!     To provide custom comments, derive a new object from CTCommentProvider that overrides
//!     the getComment function. A custom derived CTCommentProvider allows control over whether comments are 
//!     written as well as the source and formatting of the data.
//!     An easy way to create a new CTCommentProvider is by passing a function pointer, a lambda, or a "functor" object to
//!     the CTCommentProvider::new_from templated function.
class CTCommentProvider {
  public:
  	//! default constructor (does nothing in this base class)
  	CTCommentProvider(){}
	virtual ~CTCommentProvider(){}
	//! This function returns the desired comment for a given structure.
	//! The function may return "" to indicate that no comment should be 
	//! prepended to the structure's existing CT label.
	virtual std::string getComment(const structure* ct, const int structurenumber) = 0;
};

// template<typename TProvider>
// //! A class used by the 'from' function to create new templated classes that derive from CTCommentProvider.
// class ProviderTemplate : public CTCommentProvider {
//   public:
// 	const TProvider _provider;
// 	ProviderTemplate(TProvider provider) : _provider(provider) { }
// 	ProviderTemplate() { }
// 	std::string getComment(const structure*const ct, const int structurenumber) const { return _provider(ct, structurenumber); }
// };

struct CTComments {
	class NoCommentProvider : public CTCommentProvider {
		public: std::string getComment(const structure* ct, const int structurenumber) { return ""; }
	};
	class EnergyCommentProvider : public CTCommentProvider { public: std::string getComment(const structure* ct, const int structurenumber); };
	//! An implementation that never writes comments (i.e. returns "" every time)
	static NoCommentProvider None;
	//! Returns the default implementation that writes "ENERGY = GetEnergy()" if GetEnergy() != 0
	static EnergyCommentProvider Energy;
};

#endif // CT_COMMENT_PROVIDER_H