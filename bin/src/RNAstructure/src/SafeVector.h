#ifndef _SAFEVECTOR_
#define _SAFEVECTOR_

#include <cassert>
#include <vector>

//! Class derived from the STL std::vector.
template<class TYPE>
class SafeVector : public std::vector<TYPE>{
 public:

  // miscellaneous constructors
  SafeVector() : std::vector<TYPE>() {}
  SafeVector (size_t size) : std::vector<TYPE>(size) {}
  SafeVector (size_t size, const TYPE &value) : std::vector<TYPE>(size, value) {}
  SafeVector (const SafeVector &source) : std::vector<TYPE>(source) {}

};

#endif
