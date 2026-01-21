#include "arrays.h"
#include <algorithm>
#include <cfloat>
#include <iostream>
#include <type_traits>
#include "zero.h"
#include <cmath>
using std::cerr;
using std::vector;
const real_t ALMOST_ZERO = 0.0;
//initialize all the arrays to the correct size

template<typename T> Array1D<T> make_Array1D(size_t sz, T val){
	return vector<T>(sz, val);
}
template<typename T> Array2D<T> make_Array2D(size_t sz, T val){
	return vector<vector<T>>(sz,vector<T>(sz, val));
}
template<typename T> Array3D<T> make_Array3D(size_t first_dim, size_t sz, T val){
	return vector<vector<vector<T>>>(first_dim,vector<vector<T>>(sz,vector<T>(sz, val)));
}

template<typename T>
arrays<T>::arrays (int size) :
	V(make_Array3D(num_NCMs, size, IDENTITY)),
	W(make_Array2D(size, IDENTITY)),
	WL(make_Array2D(size, IDENTITY)),
	WMB(make_Array2D(size, IDENTITY)),
	WMBL(make_Array2D(size, IDENTITY)),
	W5(make_Array1D(size+1,IDENTITY)),
	W3(make_Array1D(size+1,IDENTITY)),
    inf(std::is_integral<T>::value? 10000000:0.0),
    zro(std::is_integral<T>::value? 0:1.0),
	N(size)
{
}

template<typename T>
void arrays<T>::print_V(){
	for(int p=0;p<num_pairtypes;p++)
		for(int i=0;i<N;i++)
			for(int j=0;j<N;j++)
				std::cout<<V[p][i][j]<<"\n";
}

template<typename T>
T arrays<T>::get_v(int p, int i, int j) const{
	return V[p][i][j];
}

template<typename T>
void arrays<T>::set_v(T val,int p, int i, int j){
  if(std::isnan((float) val)){
    cerr << "invalid value at "<<i<<" "<<j<<"\n";
    throw "up";
  }
	V[p][i][j] = val;
}

template<typename T>
T arrays<T>::get_w(int i, int j) const{
	return W[i][j];
}

template<typename T>
void arrays<T>::set_w(T val, int i, int j){
	W[i][j] = val;
}

template<typename T>
T arrays<T>::get_wl(int i, int j) const{
	return WL[i][j];
}

template<typename T>
void arrays<T>::set_wl(T val, int i, int j){
	WL[i][j] = val;
}

template<typename T>
T arrays<T>::get_wmb(int i, int j) const{
	return WMB[i][j];
}

template<typename T>
void arrays<T>::set_wmb(T val, int i, int j){
	WMB[i][j] = val;
}

template<typename T>
T arrays<T>::get_wmbl(int i, int j) const{
	return WMBL[i][j];
}

template<typename T>
void arrays<T>::set_wmbl(T val, int i, int j){
	WMBL[i][j] = val;
}

template<typename T>
T arrays<T>::get_w5(int i) const{
    if (i==-1) return zro;
    if(i<-1 || i>N){
        std::cerr<<"out of bounds access to w5\n";
        throw "out of bounds access";
    }
	return W5[i+1];
}

template<typename T>
void arrays<T>::set_w5(T val,int i){
    if(i<-1 || i>N){
        std::cerr<<"out of bounds access to w5\n";
        throw "out of bounds access";
    }
    if(i==-1) std::cout<<"setting w5 at i=-1\n";
	W5[i+1] = val;
}

template<typename T>
T arrays<T>::get_w3(int i) const{
    if (i==N) return zro;
    if(i<-1 || i>N){
        std::cerr<<"out of bounds access to get_w3: i="<<i<<"\n";
        throw "out of bounds access";
    }
	return W3[i];
}

template<typename T>
void arrays<T>::set_w3(T val,int i){
    if(i<-1 || i>N){
        std::cerr<<"out of bounds access to set_w3: i="<<i<<"\n";
        throw "out of bounds access";
    }
	W3[i] = val;
}

template class arrays<real_t>;
template class arrays<int>;
