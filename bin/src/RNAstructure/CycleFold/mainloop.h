#ifndef __MAINLOOP_H__
#define __MAINLOOP_H__
#include "arrays.h"
#include "sequence.h"
#include "constraints.h"
#include "NCM_parameters.h"
#include <vector>
typedef std::vector<std::vector<double>> table_t;
inline table_t make_table(size_t xdim, size_t ydim){
    return std::vector<std::vector<double>>(xdim,std::vector<double>(ydim,0.0));
}

#include "extrinsic.h"

table_t calculate_pairing_probabilities(const sequence& s, const parameters<real_t>& par, const constraints& cst, const extrinsic<real_t>& extr=extrinsic<real_t>(), const options& op=options());

void calculate_mfe(const sequence& s, const parameters<int>& par, const constraints& cst, const options& op);

template<typename T>
arrays<T> fill_arrays(const sequence& s, const parameters<T>& par,const constraints& cst, const extrinsic<T>& extr=extrinsic<T>(), const options& op=options());

template<typename T>
void calc_V_hairpin(const int d,arrays<T>& arr,const sequence& s,const parameters<T>& par,const constraints& cst);
template<typename T>
void calc_V_exterior(const int d,arrays<T>& arr,const sequence& s,const parameters<T>& par,const constraints& cst, const extrinsic<T>& extr=extrinsic<T>());
template<typename T>
void calc_V_multibranch(const int d,arrays<T>& arr,const sequence& s,const parameters<T>& par,const constraints& cst, const extrinsic<T>& extr=extrinsic<T>());
template<typename T>
void calc_V_bigloop(const int d,arrays<T>& arr,const sequence& s,const parameters<T>& par,const constraints& cst);
template<typename T>
void calc_V_stack(const int d,arrays<T>& arr,const sequence& s, const parameters<T>& par,const constraints& cst);
template<typename T>
void calc_W(const int d,arrays<T>& arr,const sequence& s, const parameters<T>& par,const constraints& cst);
template<typename T>
void calc_WL(const int d,arrays<T>& arr,const sequence& s, const parameters<T>& par,const constraints& cst);
template<typename T>
void calc_WMB(const int d,arrays<T>& arr,const sequence& s, const parameters<T>& par,const constraints& cst);
template<typename T>
void calc_WMBL(const int d,arrays<T>& arr,const sequence& s, const parameters<T>& par,const constraints& cst);
template<typename T>
void calc_W5(const int d,arrays<T>& arr,const sequence& s, const parameters<T>& par,const constraints& cst);
template<typename T>
void calc_W3(const int d,arrays<T>& arr,const sequence& s, const parameters<T>& par,const constraints& cst);
int wrap(const int j,const int N);
bool too_close_big_hairpin(const int i, const int j, const NCM_type& n);
bool too_close(const int N,const int i, const int j,const NCM_type& n, const NCM_type& m);
bool too_close(const int i, const int j, const NCM_type& h);
int ncm_length(const NCM_type& n);
bool is_double_stranded(const NCM_type& n);
void set_di_dj(const int i,const int j,int& di,int& dj, const NCM_type& n);
bool is_exterior(const int i,const int j);
bool is_interior(const int i,const int j);
template<typename T>
double coaxial(const int i,const int j,const int k, const int l, const parameters<T>& p, const sequence& s);
void show_pair_probs(const table_t& probs);
void show_pair_probs_arrayview(const table_t& probs);
template<typename T>
void show_v_array(const arrays<T>& a, const int n, const parameters<T>& p);
template<typename T>
void show_w5_w3(const arrays<T>& a, const int n);
template<typename T>
std::vector< std::pair<int,int> > probknot(const arrays<T>& a, const sequence& s,const int n, const parameters<real_t>& p);
void show_ct(const std::vector<std::pair<int,int> >& pairs, const sequence& s, const int energy=1);
template<typename T>
std::vector<std::pair<int,int>> traceback(const sequence& seq, const arrays<int>& a, const parameters<int>& p, const options& op);

//for tubofold calculation
template<typename T>
void calc_V_extrinsic(const int d, const int N, arrays<T>& arr,const extrinsic<T>& extr);
#endif
