#include "mainloop.h"
#include <iostream>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <type_traits>
#include <stack>
//#include "logdouble.h"
#include "../src/phmm/utils/xmath/log/xlog_math.h"
#include "zero.h"
using std::vector;
using std::pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::string;
const real_t Z = real_t(0.0);
const bool VERBOSE = false;

template<typename T>
bool mfe_calc() {
  return std::is_integral<T>::value;
}

bool valid_index(const int i, const int N){
  return i >= 0 && i < N;
}

bool invalid_index(const int i, const int N){
  return !valid_index(i,N);
}

template<typename T>
void show_w_arrays(const arrays<T>& arr, const int n);

double pair_probability(const arrays<real_t>& a, const int N, const sequence& s,const parameters<real_t>& p, const int i,const int j, const constraints& cst, const extrinsic<real_t>& extr=extrinsic<real_t>());

void calculate_mfe(const sequence& s, const parameters<int>& par, const constraints& cst, const options& op){
  arrays<int> arr = fill_arrays(s,par,cst,extrinsic<int>(),op);
  vector<pair<int,int>> pairs = traceback<int>(s, arr, par, op);
  if (VERBOSE){
    cout<<arr.get_w3(0)<<endl;
  }
  show_ct(pairs,s,arr.get_w3(0));
}

table_t calculate_pairing_probabilities(const sequence& s, const parameters<real_t>& par, const constraints& cst, const extrinsic<real_t>& extr, const options& op){
  const int n = s.getLength();
  table_t probs = make_table(n,n);
  arrays<real_t> pfunc = fill_arrays(s,par,cst,extr,op);
  for(int i=0;i<n;i++){
    for(int j=i+1;j<n;j++){
      probs[i][j] = probs[j][i] = pair_probability(pfunc,n,s,par,min(i,j),max(i,j),cst,extr);
    }
  }
  return probs;
}

template<typename T>
arrays<T> fill_arrays(const sequence& s,
                      const parameters<T>& par,
                      const constraints& cst,
                      const extrinsic<T>& extr,
                      const options& op/*=false*/)
{//main algorithm
  const int N = s.getLength();
  arrays<T> arr = arrays<T>(N);
  calc_W5(0,arr,s,par,cst);
  calc_W3(0,arr,s,par,cst);
  for (int d = 1; d < N; ++d){//distance between nucleotides
    calc_V_hairpin(d,arr,s,par,cst);
    calc_V_exterior(d,arr,s,par,cst,extr);
    calc_V_stack(d,arr,s,par,cst);
    calc_V_multibranch(d,arr,s,par,cst, extr);
		if(op.bigloops){
			//calc_V_bigloop(d,arr,s,par,cst);
		}
		calc_V_extrinsic(d,s.getLength(),arr,extr);
    calc_WL(d,arr,s,par,cst);
    calc_W(d,arr,s,par,cst);
    calc_WMBL(d,arr,s,par,cst);
    calc_WMB(d,arr,s,par,cst);
    calc_W5(d,arr,s,par,cst);
    calc_W3(d,arr,s,par,cst);
  }
	if(VERBOSE){
		show_v_array(arr, s.getLength(), par);
		show_w_arrays(arr,s.getLength());
    show_w5_w3(arr,s.getLength());
	}
  return arr;
}

template<typename T>
void calc_V_hairpin(const int d,arrays<T>& arr,const sequence& s,
                    const parameters<T>& par,const constraints& cst
                    ){
  if (!(d<6 && d > 1)) return;//3<=d<=6 ..
  const int N = s.getLength();
  NCM_type n;
  if(d==2) n = std::string("13");
  if(d==3) n = std::string("14");
  if(d==4) n = std::string("15");
  if(d==5) n = std::string("16");
  //find possible hairpin
  for(int i=0;i<N-d;i++){
    const int j = wrap(i+d,N);
    if (is_exterior(i,j)) continue;
		NCM_id id = NCM_id(n,s.substr(i,NCM::fivep_length(n)));
		T V;
		if(mfe_calc<T>()){
			V = conflicts(i,j,n,cst)? IDENTITY : par.energy(id);
		}
		else {
			V = conflicts(i,j,n,cst)? IDENTITY : par.energy(id);
		}
    arr.set_v(V,NCM::to_int(n),i,j);
  }
}

NCM_id coaxial_NCM_five(const int i, const int k, const std::string& s){
  char buf[] = {s.at(0),s.at(k),s.at(k+1),s.at(i)};
  std::string seq = std::string(buf);
  return NCM_id(std::string("222"),seq);
}

NCM_id coaxial_NCM_three(const int i, const int k, const std::string& s){
  char buf[] = {s.at(i),s.at(k),s.at(k+1),s.at(s.length()-1)};
  std::string seq = std::string(buf);
  return NCM_id(std::string("222"),seq);
}

template<typename T>
void calc_V_multibranch(const int d,arrays<T>& arr,const sequence& s,
                        const parameters<T>& par,const constraints& cst, const extrinsic<T>& extr)
{
  const int N = s.getLength();
  for(int i=0;i<N;i++){
    const int j = wrap(i+d,N);
    if(!allowed_pair(i,j,cst)) continue;
    auto allowed_NCMs = NCM::stack_bulge();//is_exterior(i,j)? NCM::stack_bulge() : NCM::all_NCMs();
    for(NCM_type n : allowed_NCMs){//NCM we're adding
      int di,dj;
      set_di_dj(i,j,di,dj,n);
      if(!too_close(N,i,j,n,std::string("mb"))){
        T V = arr.get_v(NCM::to_int(n),i,j);
        T W = arr.get_wmb(i+di+1,j-dj-1);
        NCM_id inner = NCM_id(i,j,n,s.toString());
				if(mfe_calc<T>()){
          V = min(V,W+par.energy(inner)+par.pair_bonus());
				}
				else {
					V += W*par.energy(inner)*par.pair_bonus()*extr.bonus(i+di,j-dj);
				}
        bool illegal = conflicts(i,j,n,cst);
        arr.set_v(illegal? IDENTITY : V,NCM::to_int(n),i,j);
      }
    }
  }
}

template<typename T>
void calc_V_bigloop(const int d,arrays<T>& arr,const sequence& s,
                    const parameters<T>& par,const constraints& cst)
				
{
  const int N = s.getLength();
  for(int i=0;i<N-d;i++){
    const int j = wrap(i+d,N);
    if(!allowed_pair(i,j,cst)) continue;
    for(NCM_type n : NCM::stack_bulge()){//NCM we're adding
			if (too_close_big_hairpin(i,j,n)) continue;
			NCM_id inner = NCM_id(i,j,n,s.toString());
      int di,dj;
      set_di_dj(i,j,di,dj,n);
			T V = arr.get_v(NCM::to_int(n),i,j);
			//big internal loop search
			for(NCM_type m: NCM::stack_bulge()){
        if(too_close(N,i,j,n,m)) continue;
				for(int ip=i+di+1;ip<j-dj-1;ip++){
					for(int jp=ip+1;jp<j-dj-1;jp++){
						T V2 = arr.get_v(NCM::to_int(m),ip,jp);
						if(mfe_calc<T>()){
							V = min(V,V2+par.energy(inner)+par.pair_bonus());
						}
						else{
							V += V2*par.energy(inner)*par.pair_bonus();
						}
					}
				}
			}
			//big hairpin loop
			if(is_interior(i,j) && !too_close_big_hairpin(i,j,n)){
				//1-NCM bigger than 6nt
				if(mfe_calc<T>()){
					V = min(V,T(0.0)+par.energy(inner)+par.pair_bonus());
				}
				else {
					V += T(1.0)*par.energy(inner)*par.pair_bonus();
				}
			}
			bool illegal = conflicts(i,j,n,cst);
			arr.set_v(illegal? IDENTITY : V,NCM::to_int(n),i,j);
    }
  }
}

template<typename T>
void calc_V_exterior(const int d,
                     arrays<T>& arr,
                     const sequence& s,
                     const parameters<T>& par,
                     const constraints& cst, const extrinsic<T>& extr)
{
  const int N = s.getLength();
  for(int i=0;i<N;i++){
    const int j = wrap(i+d,N);
    if(!is_exterior(i,j)) continue;
    for(NCM_type n : NCM::stack_bulge()){//NCM we're adding
      int di,dj;
      set_di_dj(i,j,di,dj,n);
      if(!too_close(N,i,j,n,std::string("ext"))){
        bool illegal = conflicts(i,j,n,cst);
        T V = arr.get_v(NCM::to_int(n),i,j);
				NCM_id theta = NCM_id(i,j,n,s.toString());
        if( mfe_calc<T>() ){
          T W = arr.get_w5(j-dj-1) + arr.get_w3(i+di+1);
          V = min(V,W+par.energy(theta)+par.pair_bonus());
        }
        else {
          T W = arr.get_w5(j-dj-1) * arr.get_w3(i+di+1);
          V += W*par.energy(theta)*par.pair_bonus()*extr.bonus(i+di,j-dj);
        }
        arr.set_v((illegal? IDENTITY : V),NCM::to_int(n),i,j);
      }
    }
  }
}


template<typename T>
void calc_V_stack(const int d,
                  arrays<T>& arr,
                  const sequence& s,
                  const parameters<T>& par,
                  const constraints& cst)
{
  const int N = s.getLength();
  for(int i=0;i<N;i++){
    const int j = wrap(i+d,N);
    for(NCM_type n : NCM::stack_bulge()){//NCM we're adding
      bool illegal = conflicts(i,j,n,cst);
      T V = arr.get_v(NCM::to_int(n),i,j);
      int di,dj;
      set_di_dj(i,j,di,dj,n);
      for(NCM_type m : NCM::all_NCMs()){//NCMs we're adding on to
        if(too_close(N,i,j,n,m)) continue;
				NCM_id outer = NCM_id(i,j,n,s.toString());
				NCM_id inner = NCM_id(i+di,j-dj,m,s.toString());
				T e = par.energy(outer,inner,is_exterior(i,j));
				T vp = arr.get_v(NCM::to_int(m),i+di,j-dj);
				if(conflicts(i+di,j-dj,m,cst)){
					assert(vp == IDENTITY);
				}
				if(mfe_calc<T>()){
					V = min(V,e+vp);
				}
				else {
					V += e * vp;
				}
      }
      arr.set_v(illegal? IDENTITY : V,NCM::to_int(n),i,j);
    }
  }
}

//W(i,j) = WL(i,j) + W(i,j-1)
template<typename T>
void calc_W(const int d,arrays<T>& arr,const sequence& s, const parameters<T>& par,const constraints& cst)
{
  const int N = s.getLength();
  for(int i=0;i<N;i++){
    const int j = wrap(i+d,N);
    T w = arr.get_wl(i,j);
    //check constraint
    if(allowed_unpair(j,cst)){
      if(valid_index(j-1,N)){
        if(mfe_calc<T>()){
          w  = min(w, arr.get_w(i,j-1));// * par.mb_nuc_penalty();
        }
        else{
          w += arr.get_w(i,j-1);// * par.mb_nuc_penalty();
        }
      }
    }
    arr.set_w(w,i,j);
  }
}

template<typename T>
void calc_WL(const int d,arrays<T>& arr,const sequence& s, const parameters<T>& par,const constraints& cst)
{
  const int N = s.getLength();
  for(int i=0;i<N;i++){
    const int j = wrap(i+d,N);
    T wl = IDENTITY;
    //check constraint
    if(allowed_unpair(i,cst)){
      if(valid_index(i+1,N)){
        if(mfe_calc<T>()){
          wl = min(wl, arr.get_wl(i+1,j));// * par.mb_nuc_penalty();
        }
        else {
          wl += arr.get_wl(i+1,j);// * par.mb_nuc_penalty();
        }
      }
    }
    for(NCM_type m : NCM::stack_bulge()){
      if(mfe_calc<T>()){
        wl = min(wl,arr.get_v(NCM::to_int(m),i,j));
      }
      else {
        wl += arr.get_v(NCM::to_int(m),i,j);
      }
    }
    arr.set_wl(wl,i,j);
  }
}

template<typename T>
void calc_V_extrinsic(const int d, const int N, arrays<T>& arr,const extrinsic<T>& extr)
{
  for(int i=0;i<N;i++){
    const int j = wrap(i+d,N);
    for(NCM_type n : NCM::all_NCMs()){
			T V = arr.get_v(NCM::to_int(n), i, j);
			if(mfe_calc<T>()){
				V += extr.bonus(i,j);
			}
			else {
				V *= extr.bonus(i,j);
			}
			arr.set_v(V, NCM::to_int(n), i, j);
		}
	}
}

template<typename T>
void calc_WMB(const int d,arrays<T>& arr,const sequence& s, const parameters<T>& par,const constraints& cst)
{
  const int N = s.getLength();
  for(int i=0;i<N;i++){
    //wmb(i,j)=wmb(i+1,j)*b+wmbl(i,j)
    const int j = wrap(i+d,N);//think this is wrong should be +
    T wmb = IDENTITY;
    T wmb_p = IDENTITY;
    if (allowed_unpair(j,cst) && valid_index(j-1,N))
      wmb_p = arr.get_wmb(i,j-1);
    T wmbl = arr.get_wmbl(i,j);
    if(mfe_calc<T>()){
      wmb = min(wmb_p,wmbl);
    }
    else {
      wmb = wmb_p + wmbl;
    }
    arr.set_wmb(wmb,i,j);
  }
}

template<typename T>
real_t W_coax(const int i,const int j,const arrays<T>& arr,const sequence& s,
              const parameters<T>& par){
  const int N = s.getLength();
  std::string seq = s.toString();
  real_t wcoax = 0.0;
  for(int k = i+1;wrap(k,N)<j;k++){
    for(NCM_type m : NCM::all_NCMs()){
      for(NCM_type n : NCM::all_NCMs()){
        NCM_id five(i,k,m,seq);
        NCM_id three(k,j,m,seq);
        char coaxial_nucs[] = {seq.at(i),seq.at(k),seq.at(k+1),seq.at(j)};
        NCM_id coaxial(std::string("222"),std::string(coaxial_nucs));
        real_t coax_energy = 0.0;
        wcoax += arr.get_v(NCM::to_int(n),i,k) *
          arr.get_v(NCM::to_int(m),k+1,j) *
          coax_energy;
      }
    }
  }
  return wcoax;
}

template<typename T>
void calc_WMBL(const int d,arrays<T>& arr,const sequence& s, const parameters<T>& par,const constraints& cst)
{
  const int N = s.getLength();
  for(int i=0;i<N;i++){
    const int j = wrap(i+d,N);
    T wmbl = IDENTITY;
    //check constraint
    T add_unpair = (allowed_unpair(i,cst) && valid_index(i+1,N))?arr.get_wmbl(i+1,j):IDENTITY;
    if(mfe_calc<T>()){
      wmbl = min(wmbl,add_unpair);
    }
    else {
      wmbl += add_unpair;
    }

    for(int k=i+1;wrap(k,N)!=j;k++){
      if(k==N) continue;
      int kp = wrap(k,N);
      int kpp = wrap(k+1,N);
      if(kpp==0) continue;
      //if (i==14 && j == 3){
      //cout << kp <<"\n";
        ////cout << "five"<< five<<"\n";
        ////cout << "three"<< three<<"\n";
      //}
      T five = IDENTITY;//is this the right value to start?
      for(NCM_type n : NCM::stack_bulge()){
        if(mfe_calc<T>()) {
          five = min(five,arr.get_v(NCM::to_int(n),i,kp));
        }
        else {
          five += arr.get_v(NCM::to_int(n),i,kp);
        }
      }
      if(mfe_calc<T>()) {
				T three = min(arr.get_wmbl(kpp,j), arr.get_wl(kpp,j));
        wmbl = min(wmbl, five + three);
      }
      else {
				T three = arr.get_wmbl(kpp,j)+arr.get_wl(kpp,j);
        wmbl += five * three;
      }
    }
    arr.set_wmbl(wmbl,i,j);
  }
}

template<typename T>
void calc_W5(const int d,arrays<T>& arr,const sequence& s, const parameters<T>& par,const constraints& cst)
{
  T w5 = IDENTITY;
  if(allowed_unpair(d,cst)){
    T add_unpair = arr.get_w5(d-1);
    if(mfe_calc<T>()) {
      w5 = min(w5, add_unpair);
    }
    else {
      w5 += add_unpair;
    }
  }
  for(NCM_type n:NCM::stack_bulge()){
    T v = arr.get_v(NCM::to_int(n),0,d);
    if(mfe_calc<T>()) {
      w5 = min(w5, v);
    }
    else {
      w5 += v;
    }
  }
  for(int k=0;k<d;k++){
    for(NCM_type n:NCM::stack_bulge()){
      T w5_p = arr.get_w5(k);
      T v = arr.get_v(NCM::to_int(n),k+1,d);
      if(mfe_calc<T>()) {
        w5 = min(w5, w5_p + v);
      }
      else {
        w5 += w5_p*v;
      }

    }
  }
  arr.set_w5(w5,d);
}

template<typename T>
void calc_W3(const int d,arrays<T>& arr,const sequence& s, const parameters<T>& par,const constraints& cst)
{
  const int N = s.getLength();
  T w3 = IDENTITY;
  if(allowed_unpair(N-d-1,cst)){
    T add_unpaired = arr.get_w3(N-d);
    if(mfe_calc<T>()){
      w3 = min(w3,add_unpaired);
    }
    else {
      w3 += add_unpaired;
    }
  }
  for(NCM_type n:NCM::stack_bulge()){
    T v  = arr.get_v(NCM::to_int(n),N-d-1,N-1);
    if(mfe_calc<T>()){
      w3 = min(w3, v);
    }
    else {
      w3 += v;
    }
  }
  for (int k=N-1;k>N-d-1;k--){
    for(NCM_type n:NCM::stack_bulge()){
      T v = arr.get_v(NCM::to_int(n),N-d-1,k-1);
      T w3_p = arr.get_w3(k);
      if(mfe_calc<T>()){
        w3 = min(w3, v+w3_p);
      }
      else {
        w3 += v*w3_p;
      }
    }
  }
  arr.set_w3(w3,N-d-1);
}

int wrap(const int j,const int N)
{
  assert(j<2*N);
  return j<N? j : j-N;
}

bool is_interior(const int i, const int j)
{
  assert(i!=j);
  return i<j;
}

bool is_exterior(const int i, const int j)
{
  assert(i!=j);
  return i>j;
}


bool too_close_big_hairpin(const int i, const int j, const NCM_type& n){
	int di, dj;
	set_di_dj(i,j,di,dj,n);
	return (j-dj)-(i+di)+1 <= 6;
}
//double strand NCM
//how to use:
//if exterior: n is the NCM being added on and it's on the OUTSIDE
//if interior: n is the NCM bring added on and it's on the INSIDE
bool too_close(const int N,const int i, const int j,
               const NCM_type& n, const NCM_type& m)
{
  if (i<0 || i>N-1 || j<0 || j>N-1) return true;
  const int HINGE_LENGTH = 2;
  if(m==std::string("mb")){
    if(!is_double_stranded(n)){
      return i==0 || j==N-1 || j-i != ncm_length(n)-1;
    }
    else{
			int di,dj;
			set_di_dj(i,j,di,dj,n);
      if(is_interior(i,j)){
        return (j-i+1) < ncm_length(n)+2;//+2 so that there's a wmb fragment inbetween
      }
      else{
        return i+di > N-2 || j-dj < 1;//have to have nucs on the end
      }
    }
    return false;
  }
  if(m==std::string("big")){
    //TODO add a test for this case
    assert(is_double_stranded(n));
		int di,dj;
		set_di_dj(i,j,di,dj,n);
		if(is_interior(i,j)){
			return (j-i+1) < ncm_length(n)+2;
      //+2 so that there's a fragment inbetween
      //should consider the possibility of adding 4 more for smallest possible stem
      //(212 and 13)+1 for an unpaired nuc
		}
		else{
      //			return (i+di-N-1) + (j-dj) 
			return i+di > N-2 || j-dj < 1;//have to have nucs on the end
		}
    return false;
  }
  if(m==std::string("ext")){
		if (is_interior(i,j)) throw "index error\n";
    const int dj = NCM::fivep_length(n)-1;
    const int di = NCM::threep_length(n)-1;
    return (j-dj<0) || (i+di>N-1);
  }
  if(is_interior(i,j)){
    const int dij = j-i+1;//std::max(i,j)-std::min(i,j); 6
    const int dnm = ncm_length(m)+ncm_length(n) - HINGE_LENGTH;//we have to subtract 2 for the hinge 6 
    if(!is_double_stranded(n) || !is_double_stranded(m))
      return dij != dnm;
    return dij<dnm;
  }
  else{//if is_exterior(i,j)
    const int di = NCM::threep_length(m)+NCM::threep_length(n)-2;//-2 becase i and j are in the hinge
    const int dj = NCM::fivep_length(m)+NCM::fivep_length(n)-2;
    return (j-dj<0) ||  (i+di>N-1);
  }
}

//single strand NCM
bool too_close(const int i, const int j, const NCM_type& h){
  const int dij = std::max(i,j)-std::min(i,j);
  return dij<ncm_length(h);
}

int ncm_length(const NCM_type& n){
  return NCM::fivep_length(n)
    + (is_double_stranded(n)? NCM::threep_length(n):0);
}

bool is_double_stranded(const NCM_type& n){
  return n[0]-'0'==2;
}

void set_di_dj(const int i,const int j,int& di,int& dj, const NCM_type& n){
  const int fp = NCM::fivep_length(n) - 1;//adjust by 1 for nuc in the hinge
  const int tp = NCM::threep_length(n) - 1;
  di = is_interior(i,j)?fp:tp;
  dj = is_interior(i,j)?tp:fp;
}

double pair_probability(const arrays<real_t>& a, const int N, const sequence& s,const parameters<real_t>& p, const int i, const int j, const constraints& cst, const extrinsic<real_t>& extr){
  double prob = 0.0;
  if (i==j) return prob;
  const real_t Q = a.get_w3(0);
  //probability that the pair is in a hinge
  for(NCM_type n : NCM::all_NCMs()){ //interior NCM
    for(NCM_type m : NCM::stack_bulge()){ //exterior NCM
      const int ip = i + NCM::fivep_length(n)-1;
      const int jp = j - NCM::threep_length(n)+1;
      const int ipp = i - NCM::fivep_length(m)+1;//+1 for the hinge
      const int jpp = j + NCM::threep_length(m)-1;
      if(ip> N-1 || jp < 0) continue;
      if(!is_double_stranded(n))
        if(!(j-i+1==ncm_length(n)))
          continue;

      if(ipp<0 || jpp > N-1 || too_close(N,ipp,jpp,n,m)){
        continue;
      }
      if(!is_double_stranded(n)){
      }
      const real_t v = a.get_v(NCM::to_int(n),i,j);
      const real_t vp = a.get_v(NCM::to_int(m),j,i);
      if(conflicts(i,j,n,cst) || conflicts(j,i,m,cst)) continue;
      if(v==real_t(0.0) || vp == real_t(0.0)) continue;
      NCM_id inner = NCM_id(i,j,n,s.toString());
      NCM_id outer = NCM_id(ipp,jpp,m,s.toString());
      //real_t K = 0.0;//p.junction_energy(outer,inner,false);
      real_t K = p.junction_energy(outer,inner,false) / (extr.bonus(i,j));
      /*
      real_t a = p.junction_energy(outer,inner,false);
      real_t b = p.junction_energy(inner,outer,true);
      real_t d = max(a,b) - min(a,b);
      cout <<"equal? "<< (d < real_t(0.00001));
      */
      prob += (double) ((v*vp*K) / Q);
      if (v*vp*K/Q > 1.0){
        //cout<< "PROB ERROR n"<<n<<" m "<<m<<" i "<<i<<" j "<<j<<" v "<<v<<" vp "<<vp<<" K "<<K<<" Q "<<Q<<" prob "<<(double)((v*vp*K)/Q)<<" "<<endl;

      }
      if (prob>1.0){
        //cout<<"PROBABILITY ERROR "<<n<<" "<<m<<" "<<i<<" "<<j<<" "<<v<<" "<<vp<<" "<<K<<" "<<Q<<" "<<(double)((v*vp*K)/Q)<<" "<<prob<<endl;
      }
    }
  }
  //possibility that this is an exterior loop
  bool ext = true;
  bool mb = true;
	bool big_loop = false;
  if(ext){
    for(NCM_type n : NCM::stack_bulge()){//NCM we're adding
      const int ip = std::min(i,j);
      const int jp = std::max(i,j);
      if(too_close(ip,jp,n)) continue;
      NCM_id inner = NCM_id(ip,jp,n,s.toString());
      real_t v = a.get_v(NCM::to_int(n),ip,jp);
      real_t vp = a.get_w5(ip-1) * a.get_w3(jp+1);
      prob += (double) ((v*vp) / Q);
    }
  }
  //possibility that this pair closes a multibranch loop
  if(mb){
    int ip = std::min(i,j);
    int jp = std::max(i,j);
    for(NCM_type n : NCM::stack_bulge()){//NCM we're adding
      int ipp = ip + NCM::fivep_length(n)-1;
      int jpp = jp - NCM::threep_length(n)+1;
      if(ipp>=jpp) continue;
      //mb loop is on the interior
      if (ip>0 && jp < N-1 && !too_close(N,ipp,jpp,n,"mb")){
        real_t v = a.get_wmb(ip+1,jp-1); 
        real_t vp = a.get_v(NCM::to_int(n),jp,ip);
        prob += (double) ((v*vp) / Q);
      }
      //mb loop is on the exterior
      if(!too_close(N,jpp,ipp,n,"mb")){
        real_t v = a.get_v(NCM::to_int(n),ip,jp);
        real_t vp = a.get_wmb(jp+1,ip-1);
        prob += (double) ((v*vp) / Q);
      }
    }
  }
	if(big_loop){
    int ip = std::min(i,j);
    int jp = std::max(i,j);
    for(NCM_type n : NCM::stack_bulge()){//NCM we're adding
			if(!too_close(N,ip,jp,n,"big")){
				//iloop
        real_t v = a.get_w(ip-1,jp+1);
        real_t vp = a.get_v(NCM::to_int(n),jp,ip);
        prob += (double) ((v*vp) / Q);
				//hairpin
				v = 1.0;
        prob += (double) ((v*vp) / Q);
			}
			if(!too_close(N,jp,ip,n,"big")){
        real_t v = a.get_v(NCM::to_int(n),ip,jp);
        real_t vp = a.get_w(jp-1,ip+1);
        prob += (double) ((v*vp) / Q);
			}
		}
	}
	//assert(prob <= 1.0000000001);
	//if(VERBOSE) cout<<"final prob "<<prob<<endl;
  return (double) prob;
}

void show_pair_probs(const table_t& probs){
  cout <<"probabilities"<<endl;
  for(size_t i=0;i<probs.size();i++){
    for(size_t j=0;j<probs.size();j++){
      cout<<i<<"\t"<<j<<"\t"<<probs[i][j]<<endl;
    }
  }
  cout<<endl;
}

void show_pair_probs_arrayview(const table_t& probs){
  int start_index = 1; // change to 0 for 0-indexed arrays
  cout <<"probs"<<endl;
  for(size_t i=0;i<probs.size();i++){
    cout<<"\t"<<(i+start_index);
  }
  cout<<endl;
  for(size_t i=0;i<probs.size();i++){
    cout<<(i+start_index);
    for(size_t j=0;j<probs.size();j++){
      cout<<"\t"<<probs[i][j];
    }
    cout<<endl;
  }
  cout<<endl;
}

template<typename T>
vector< pair<int,int> > probknot(const arrays<T>& a, const sequence& s,const int n, const parameters<real_t>& p){
  vector< vector <double> > t(n,vector<double>(n)); //pair probabilities table
  //fill with pair probs
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      t[i][j] = pair_probability(a,n,s,p,std::min(i,j),std::max(i,j));
    }
  }
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      assert(t[i][j] - t[j][i] < 1e-08);
    }
  }
  vector<double> best(n, 0.0);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if(t[i][j] > best[i]) best[i] = t[i][j];
    }
  }
  vector<int> pr;
  for(int i=0;i<n;i++){
    pr.push_back(i);
    for(int j=0;j<n;j++){
      if(t[i][j] > t[i][pr[i]]){
        pr[i] = j;
      }
    }
  }
  std::vector< std::pair<int,int> > v;
  for(int i=0;i<n;i++){
    if (pr[i] != i && i<pr[i] && pr[pr[i]]==i)
      v.push_back(std::make_pair(i,pr[i]));
  }
  return v;
}

void show_ct(const vector<pair<int,int>>& pairs, const sequence& s, const int energy){
	string e = energy==1? string("") : 
		string("\tenergy: ")
		+std::to_string(energy / 10000)
		+ string(".")
		+std::to_string(-energy%10000);
  //get the pairing info
  vector<int> v(s.getLength());
  for(int i=0;i<s.getLength();i++){
    v[i] = i;
  }
  for(const pair<int,int> p : pairs){
    v[p.first] = p.second;
    v[p.second] = p.first;
  }

  cout<< s.getLength()<<"\t"<<s.getTag()<<e<<endl;
  for(int i=0;i<s.getLength();i++){
    const int j = v[i];
    cout<< i+1<< "\t"
        << s.toString()[i] <<"\t"
        << i << "\t"
        << i+2 << "\t"
        <<(j!=i? j+1 : 0) << "\t"
        << i+1
        << endl;
  }
}

typedef struct traceback_pair{
  int i;
  int j;
  int energy;
  bool paired;
} traceback_t;

traceback_t mkPair(const int i, const int j, const int energy, const bool paired, bool verbose=VERBOSE){
	assert(i<j);
	if(verbose) cout<<"pushed a pair "<<i<<" "<<j<<" "<<energy<<endl;
  traceback_t t;
  t.i = i;
  t.j = j;
  t.energy = energy;
  t.paired = paired;
  return t;
}

void print_stack(const std::stack<traceback_t>& stk){
  std::stack<traceback_t> cpy = stk;
  cout<<"\tstate of stack:\n";
  while(!cpy.empty()){
    traceback_t tt = cpy.top();
    cout<<"\t"<<tt.energy<<" "<<tt.i<<" "<<tt.j<<endl;
    cpy.pop();
  }

}

template<typename T>
vector<pair<int,int>> traceback(const sequence& seq, const arrays<int>& a, const parameters<int>& par, const options& op){
	bool verbose=VERBOSE;
  const int n = a.N;
  vector<pair<int,int>> pairs;
	if(a.get_w3(0) == 0){
		return pairs;
	}
  std::stack<traceback_t> s;
  s.push(mkPair(0,n-1,a.get_w3(0),false));
  if (verbose) { cout << "energy "<<a.get_w3(0)<<endl;}
  while (!s.empty()){
		//print_stack(s);
    const traceback_t t = s.top();
    s.pop();
    //nucs close a hairpin
		if(true) {
			for(NCM_type h : NCM::hairpins()){
				//check distance between nucs is the size of the hairpin
				if(t.j - t.i + 1 == ncm_length(h)){
					const T v_hairpin = t.energy;
					NCM_id hp = NCM_id(h, seq.substr(t.i,ncm_length(h)));
					if(v_hairpin == par.energy(hp)){
						//found a hairpin loop
						//put the pair in
						//nothing new on the stack
						pairs.push_back({t.i, t.j});
						if (verbose) cout<<"found a hairpin\n";
						if(verbose) {
							T e = par.energy(hp,true);
							cout<<"total energy: "<<e<<endl;
						}
						goto found;
					}
				}
			}
			//nucs close an NCM
			for(NCM_type n : NCM::stack_bulge()){
				//const T v_stack = a.get_v(NCM::to_int(n), t.i, t.j);
				const T v_stack = t.energy;
				const int ip = t.i + NCM::fivep_length(n) - 1;
				const int jp = t.j - NCM::threep_length(n) + 1;
				for (NCM_type m : NCM::all_NCMs()){
					if(too_close(seq.getLength(),t.i,t.j,n,m)) continue;
          //TODO can I just check V like this? what if value is 0?
					const T vp = a.get_v(NCM::to_int(m),ip,jp);
					NCM_id stack = NCM_id(t.i, t.j, n, seq.toString());
					NCM_id mid = NCM_id(ip, jp, m, seq.toString());
					//check that t.energy is the energy for stacking
					//if so, add this pair and push the next pair on the stack
          //can I add in the pairs for m, too??
          //cout << par.energy(stack, m, p, false) + vp<<endl;
          if(v_stack == par.energy(stack,mid,false) + vp){
            if (verbose) cout<<"found a stack\n";
            if(verbose) {
              T e = par.energy(stack,mid,false,true);
              cout<<"total energy: "<<e<<endl;
            }
							
            pairs.push_back({t.i,t.j});
            s.push(mkPair(ip,jp,vp,true));
            goto found;
					}
				}
        //stacking on a multibranch loop
        if(!too_close(seq.getLength(), t.i, t.j, n, "mb")){
          NCM_id on_mb = NCM_id(t.i,t.j,n,seq.toString());
          T mb = a.get_wmb(ip+1,jp-1)+par.energy(on_mb)+par.pair_bonus();
          if(v_stack /*== a.get_wmb(ip,jp)*/ == mb){
            if (verbose) cout<<"found a stack on a multibranch\n";
            pairs.push_back({t.i,t.j});
            pairs.push_back({ip,jp});
            //s.push(mkPair(ip,jp,a.get_wmb(ip,jp),true));
            s.push(mkPair(ip+1, jp-1, a.get_wmb(ip+1,jp-1),false));
            goto found;
          }
        }
        if(op.bigloops){
					//stacking on a big internal loop
          if(!too_close(seq.getLength(), t.i, t.j, n, "mb")){
            NCM_id on_bl = NCM_id(t.i,t.j,n,seq.toString());
            T bl = a.get_w(ip+1,jp-1)
              +par.energy(on_bl)+par.pair_bonus();
            if(v_stack == bl){
              if (verbose) cout<<"found a big internal loop\n";
              pairs.push_back({t.i,t.j});
              pairs.push_back({ip,jp});
              s.push(mkPair(ip+1, jp-1, a.get_w(ip+1,jp-1),false));
              goto found;
            }
          }
        }
				if(op.bigloops){
					//stacking on a big hairpin loop
          if(!too_close_big_hairpin(t.i,t.j,n)){
            NCM_id on_bl = NCM_id(t.i,t.j,n,seq.toString());
            T bl = 0 + par.energy(on_bl)+par.pair_bonus();
            if(v_stack == bl){
              if (verbose) cout<<"found a big hairpin at "<<t.i<<" "<<t.j<<endl;
              pairs.push_back({t.i,t.j});
              pairs.push_back({ip,jp});
              goto found;
            }
					}
				}
			}
		} //end if t.paired (currently if true)

		//bifurcation
		//if i closed a stem with another nucleotide between i and j,
		//push the stem and the rest of the multibranch loop
		for(int k=t.i+1;k<t.j;k++){
			//case with more than two stems
			if(t.energy == a.get_wl(t.i,k)+a.get_wmb(k+1,t.j)){
				//pairs.push_back({t.i,k});
				s.push(mkPair(k+1, t.j, a.get_wmb(k+1,t.j),false));
				s.push(mkPair(t.i, k, a.get_wl(t.i,k),false));
				if (verbose) cout<<"found a bifurcation\n";
				goto found;
			}
			//case with exactly two stems
			if(t.energy == a.get_wl(t.i,k)+a.get_w(k+1,t.j)){
				//pairs.push_back({t.i,k});
				s.push(mkPair(k+1, t.j, a.get_w(k+1,t.j),false));
				s.push(mkPair(t.i, k, a.get_wl(t.i,k),false));
				if (verbose) cout<<"found a bifurcation\n";
				goto found;
			}
		}
    //one of the nucs is unpaired in a multibranch loop

		//if this energy is the same as an adjacent value of w or wmb
		//add no pair, push that pair "1 in"

		//the nuc is unpaired on the 5' side
		//single stem case
		if(t.energy == a.get_w(t.i+1,t.j)){
      if (verbose) cout<<"found an unpaired nuc 5'\n";
			s.push(mkPair(t.i+1, t.j, a.get_w(t.i+1,t.j),false));
			goto found;
		}
		//multiple stem case
		if(t.energy == a.get_wmb(t.i+1,t.j)){
      if (verbose) cout<<"found an unpaired nuc 5'\n";
			s.push(mkPair(t.i+1, t.j, a.get_wmb(t.i+1,t.j),false));
			goto found;
		}
		//the nuc is unpaired on the 3' side
		//single stem case
		if(t.energy == a.get_w(t.i,t.j-1)){
      if (verbose) cout<<"found an unpaired nuc 3'\n";
			s.push(mkPair(t.i, t.j-1, a.get_w(t.i,t.j-1),false));
			goto found;
		}
		//multiple stem case
		if(t.energy == a.get_wmb(t.i,t.j-1)){
      if (verbose) cout<<"found an unpaired nuc 3'\n";
			s.push(mkPair(t.i, t.j-1, a.get_wmb(t.i,t.j-1),false));
			goto found;
		}
		throw "traceback failed\n";
  found: ;
		//std::cout <<"found a pair between "<<pairs.back().first<< " "<<pairs.back().second<<"\n";
  }
  return pairs;
}

template<typename T>
void show_v_array(const arrays<T>& a, const int n, const parameters<T>& p){

  std::cout<<std::setprecision(4);
  for(NCM_type m : NCM::all_NCMs()){
    std::cout<<m<<" ARRAY"<<std::endl<<"\t ";
    for(int i=0;i<n;i++){ std::cout<<i<<"\t";}
    std::cout<<std::endl;
    for(int i=0;i<n;i++){
      std::cout<<i<<"\t";
      for(int j=0;j<n;j++){
        std::cout<<" "<<a.get_v(NCM::to_int(m),i,j)<<"\t";
      }
      std::cout<<std::endl;
    }
  }
}

template<typename T>
void show_w_arrays(const arrays<T>& a, const int n){
  std::cout<<std::setprecision(4);
	std::cout<<"w ARRAY"<<endl<<"\t ";
	for(int i=0;i<n;i++){ std::cout<<i<<"\t";}
	std::cout<<std::endl;
	for(int i=0;i<n;i++){
		std::cout<<i<<"\t";
		for(int j=0;j<n;j++){
      std::cout<<" "<<a.get_w(i,j)<<"\t";
    }
		std::cout<<std::endl;
	}
	std::cout<<"wl ARRAY"<<endl<<"\t ";
	for(int i=0;i<n;i++){ std::cout<<i<<"\t";}
	std::cout<<std::endl;
	for(int i=0;i<n;i++){
		std::cout<<i<<"\t";
		for(int j=0;j<n;j++){
      std::cout<<" "<<a.get_wl(i,j)<<"\t";
    }
		std::cout<<std::endl;
	}
	std::cout<<"wmb ARRAY"<<endl<<"\t ";
	for(int i=0;i<n;i++){ std::cout<<i<<"\t";}
	std::cout<<std::endl;
	for(int i=0;i<n;i++){
		std::cout<<i<<"\t";
		for(int j=0;j<n;j++){
      std::cout<<" "<<a.get_wmb(i,j)<<"\t";
    }
		std::cout<<std::endl;
	}
	std::cout<<"wmbl ARRAY"<<endl<<"\t ";
	for(int i=0;i<n;i++){ std::cout<<i<<"\t";}
	std::cout<<std::endl;
	for(int i=0;i<n;i++){
		std::cout<<i<<"\t";
		for(int j=0;j<n;j++){
      std::cout<<" "<<a.get_wmbl(i,j)<<"\t";
    }
		std::cout<<std::endl;
	}
}

template<typename T>
void show_w5_w3(const arrays<T>& a, const int n){
  std::cout<<"w5 ARRAY"<<std::endl;
  for(int i=0;i<n;i++)
    std::cout<<i<<"\t";
  std::cout<<std::endl;
  for(int i=0;i<n;i++)
    std::cout<<a.get_w5(i)<<"\t";
  std::cout<<std::endl;
  std::cout<<"w3 ARRAY"<<std::endl;
  for(int i=0;i<n;i++)
    std::cout<<i<<"\t";
  std::cout<<std::endl;
  for(int i=0;i<n;i++)
    std::cout<<a.get_w3(i)<<"\t";
  std::cout<<std::endl;

}
