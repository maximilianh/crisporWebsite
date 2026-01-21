/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 *               Michael Sloma, 2015
 */

#include <iostream>
#include "DynProgArray.h"
#include "defines.h"

#ifdef EXTENDED_DOUBLE
#include "extended_double.h"
#endif

// this "function" checks at compile time whether
// two types are the same (i.e. one is an alias
// of the other)
template<typename T, typename U>
struct is_same {
	enum { value = 0 };
};

template<typename T>
struct is_same<T, T> {
	enum { value = 1 };
};

// this chooses a reasonable default value for "infinite
// energy". T should be either PFPRECISION or integersize.
template<typename T>
int get_infinity () {
	if (is_same<T,PFPRECISION>::value){
		return ZERO;
	}
	else if (is_same<T,integersize>::value){
		return INFINITE_ENERGY;
	}
	else {
		std::cerr<<"warning: infinite energy not set in DynProgArray for this type\n";
		return INFINITE_ENERGY;
	}
}

template<typename T>
DynProgArray<T>::DynProgArray(int size, int inf /* = -1 */) {
	if(inf==-1){
		infinite = (T) get_infinity<T>();
	}
	else {
		infinite = (T) inf;
	}

	Size = size;
	int i,j;
	dg = new T *[size+1];

	for (i=0;i<=(size);i++){
		dg[i] = new T [size+1];
	}
	for (i=0;i<=size;i++){
		for (j=0;j<size+1;j++){
			dg[i][j] = infinite;
		}
	}

	//Now move pointers, to facilitate fast access by
    //avoiding arithmetic during array access function
	for (i=0;i<=size;++i){
		dg[i]-=i;
	}
    //columns are i index, rows are j index
    //because i>j, array is now shaped like this:
    //     n
    // |    |
    //  |    |
    //   |    |
    //    |    |
    //         2n
    //j>n means this is an exterior fragment

}

// the destructor deallocates the space used
template<typename T>
DynProgArray<T>::~DynProgArray(){
	for (int i=0;i<=Size;i++){
		//move pointers back before deleting
		dg[i]+=i;

		//now delete
		delete[] dg[i];
	}
	delete[] dg;
}

template<typename T>
DynProgArrayT<T>::DynProgArrayT(int size, int inf) {
	if(inf==-1){
		infinite = (T) get_infinity<T>();
	}
	else {
		infinite = (T) inf;
	}
	Size = size;
	dg = new T *[2*size+1];

    //because this is transpose of DynProgArray,
    //i index is rows, j index is columns
    //array is shaped like this:
    // ||
    // | |
    // |  |
    // |   |  n
    //  |   |
    //   |  |
    //    | |
    //     || 2n
    //      n
	for (int i=0;i<=2*size;i++){
        int rowlength = i<=size ? i+1 : 2*size+1-i;
		dg[i] = new T [rowlength];
        for (int j=0;j<rowlength;j++){
            dg[i][j] = infinite;
        }
	}
    //move pointer so we don't have to do arithmetic during array access
    for(int i = size+1;i<=2*size;i++){
        dg[i] -= i-size;
    }
}

// the destructor deallocates the space used
template<typename T>
DynProgArrayT<T>::~DynProgArrayT(){
	for (int i=0;i<=2*Size;i++){
        if (i>Size){
            dg[i] += (i-Size);
        }
		delete[] dg[i];
	}
	delete[] dg;
}


template<typename T>
DynProgArrayU<T>::DynProgArrayU(int size, int inf /* = -1 */) {
	
	if(inf==-1){
		infinite = (T) get_infinity<T>();
	}
	else {
		infinite = (T) inf;
	}

	Size = size;
	int i,j;
	dg = new T *[size];

	for (i=0;i<(size);i++){
		dg[i] = new T [size-i];
	}

	for (i=0;i<size;i++){
		for (j=0;j<size-i;j++){
			dg[i][j] = infinite;
		}
	}

	//Now move pointers, to facilitate fast access by
    //avoiding arithmetic during array access function
	for (i=0;i<size;++i){
		dg[i]-=i;
	}
    //columns are i index, rows are j index
    //because i>j, array is now shaped like this:
    //     n
    // |   |
    //  |  |
    //   | |
    //    ||

}

// the destructor deallocates the space used
template<typename T>
DynProgArrayU<T>::~DynProgArrayU(){
	for (int i=0;i<Size;i++){
		//move pointers back before deleting
		dg[i]+=i;

		//now delete
		delete[] dg[i];
	}
	delete[] dg;
}


template class DynProgArray<int>;
template class DynProgArray<long>;
template class DynProgArray<long long>;
template class DynProgArray<short>;
template class DynProgArray<float>;
template class DynProgArray<double>;
template class DynProgArray<long double>;
template class DynProgArray<log_double>;
template class DynProgArrayT<int>;
template class DynProgArrayT<long>;
template class DynProgArrayT<long long>;
template class DynProgArrayT<short>;
template class DynProgArrayT<float>;
template class DynProgArrayT<double>;
template class DynProgArrayT<long double>;
template class DynProgArrayT<log_double>;
template class DynProgArrayU<int>;
template class DynProgArrayU<long>;
template class DynProgArrayU<long long>;
template class DynProgArrayU<short>;
template class DynProgArrayU<float>;
template class DynProgArrayU<double>;
template class DynProgArrayU<long double>;
template class DynProgArrayU<log_double>;


#ifdef EXTENDED_DOUBLE
template class DynProgArray<extended_double>;
template class DynProgArrayT<extended_double>;
template class DynProgArrayU<extended_double>;

#endif
