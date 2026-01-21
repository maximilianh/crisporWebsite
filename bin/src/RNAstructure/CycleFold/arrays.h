#ifndef __ARRAYH__
#define __ARRAYH__
#include "constants.h"
#include <vector>

template<typename T> using Array1D = std::vector<T>;
template<typename T> using Array2D = std::vector<std::vector<T>>;
template<typename T> using Array3D = std::vector<std::vector<std::vector<T>>>;

template<typename T>
class arrays{
  private:
  	//each array, first dimension is interior/exterior
  	//last two are indices in RNA
  	//for V, second dimension is pair type
  	//square arrays real_ts memory consumption but it's only a few
  	//MB at most so whatever
	Array3D<T> V;
	Array2D<T> W;
	Array2D<T> WL;
	Array2D<T> WMB;
	Array2D<T> WMBL;
	Array1D<T> W5;
	Array1D<T> W3;
    const T inf;
    const T zro;

  public:
	arrays(int size);
  	const int N;
	T get_v(int p,int i,int j) const;
	void set_v(T val,int p,int i,int j);
	T get_w(int i,int j) const;
	void set_w(T val,int i,int j);
	T get_wl(int i,int j) const;
	void set_wl(T val,int i,int j);
	T get_wmb(int i,int j) const;
	void set_wmb(T val,int i,int j);
	T get_wmbl(int i,int j) const;
	void set_wmbl(T val,int i,int j);
	T get_w5(int i) const;
	void set_w5(T val,int i);
	T get_w3(int j) const;
	void set_w3(T val,int j);
	void print_V();
};

template<typename T> void fill_array(Array1D<T>& arr, T);
template<typename T> void fill_array(Array2D<T>& arr, T);
template<typename T> void fill_array(Array3D<T>& arr, T);

#endif
