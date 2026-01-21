#pragma once

#ifndef ARRAYCLASS_H
#define ARRAYCLASS_H

#include <cassert>
#include "defines.h"
#include "pfunction_math.h"

// DynProgArray ("dynamic programming array")
// encapsulates the large 2-d arrays  used by
// the various dynamic programming algorithms

template <typename T>
class DynProgArray {
  public:
	int Size;
	T **dg;
	T infinite;

	// the constructor allocates the space needed by the arrays
	// the inf parameter specifies what the array uses to represent
	// "infinite energy".
	//
	// the constructor will choose a sensible default value if the template is
	// instantiated for the same type as PFPRECISION or integersize (default
	// double and short int, respectively)
	//
	// if it is instantiated for other types, a reasonable default must
	// be provided to avoid the possibility of surprising behavior.
	DynProgArray(int size, int inf = -1);

	// the destructor deallocates the space used
	~DynProgArray();

	// f is an integer function that references the correct element of
	// the array
	T &f(int i, int j);
	const T &f(int i, int j) const;
};

template <typename T>
inline T& DynProgArray<T>::f(int i, int j) {
	if (i > j) {
		return infinite;
	}
	if (i > Size) {
		i -= Size;
		j -= Size;
	}
	return dg[i][j];
}

template <typename T>
inline const T& DynProgArray<T>::f(int i, int j) const {
	if (i > j) {
			return infinite;
	}
	if (i > Size) {
			i -= Size;
			j -= Size;
	}
	return dg[i][j];
}

// DynProgArrayT has an identical interface to DynProgArray
// but its layout in memory has been transposed
// this allows for efficient traversal of columns
template<typename T>
class DynProgArrayT {
  public:
	int Size;
	T **dg;
	T infinite;

	// the constructor allocates the space needed by the arrays
	DynProgArrayT(int size, int inf = -1);

	// the destructor deallocates the space used
	~DynProgArrayT();

	// f is an integer function that references the correct element of
	// the array
	T& f(int i, int j);
	const T& f(int i, int j) const;
};

template<typename T>
inline T& DynProgArrayT<T>::f(int i, int j) {
	if (i > Size) {
		i -= Size;
		j -= Size;
	}
	if (i > j) {
		return infinite;
	}
	return dg[j][i];
}

template <typename T>
inline const T& DynProgArrayT<T>::f(int i, int j) const {
	if (i > Size) {
			i -= Size;
			j -= Size;
	}
	if (i > j) {
			return infinite;
	}
	return dg[j][i];
}

template<typename F, typename T>
inline void copyDPArray(const F& from, T& to){
	assert(from.Size==to.Size);
	for(int i=1;i<from.Size;i++){
		for(int j=i+1;j<=i+from.Size;j++){
			to.f(i,j) = from.f(i,j);
		}
	}
}
#endif

template <typename T>
class DynProgArrayU {
  public:
	int Size;
	T **dg;
	T infinite;

	// the constructor allocates the space needed by the arrays
	// the inf parameter specifies what the array uses to represent
	// "infinite energy".
	//
	// the constructor will choose a sensible default value if the template is
	// instantiated for the same type as PFPRECISION or integersize (default
	// double and short int, respectively)
	//
	// if it is instantiated for other types, a reasonable default must
	// be provided to avoid the possibility of surprising behavior.
	DynProgArrayU(int size = 1, int inf = -1);

	// the destructor deallocates the space used
	~DynProgArrayU();

	// f is an integer function that references the correct element of
	// the array
	T &f(int i, int j);
	const T &f(int i, int j) const;
};

template <typename T>
inline T& DynProgArrayU<T>::f(int i, int j) {
	if (i > j) {
		return infinite;
	}
	if (i > Size) {
		i -= Size;
		j -= Size;
	}
	return dg[i][j];
}

template <typename T>
inline const T& DynProgArrayU<T>::f(int i, int j) const {
	const T& ret = this->f(i,j);
	return ret;
}
