#include "extrinsic.h"
#include "mainloop.h"
#include "zero.h"
#include <cassert>
//#include "logdouble.h"
#include "../src/phmm/utils/xmath/log/xlog_math.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;

template<typename T>
extrinsic<T>::extrinsic() : initialized(false), gamma(0.3) {}

template<typename T>
extrinsic<T>::extrinsic(const sequence& s,
                     const vector<const sequence*> seqs,
                     const vector<const table_t*> probs,
                     const vector<const alignment*> alignments,
					 const T gamma)
: initialized(true),
  gamma(gamma)
{
  //cout<<"CALCULATING EXTRINSIC INFO\n";
  //for (auto s : seqs) cout << s->toString()<<"\n";
  //for (auto p : probs) show_pair_probs_arrayview(*p);
  //for (auto a : alignments) a->show();
  //cout<<"gamma "<<gamma<<"\n";
	bonuses = vector<vector<T>>(s.getLength(),vector<T>(s.getLength(),0.0));
  assert((seqs.size() == alignments.size()) && (alignments.size() == probs.size()));
    auto seq = seqs.begin();
    auto aln = alignments.begin();
    auto prob = probs.begin();
    for( ; seq!=seqs.end() ; ++seq,++aln,++prob){
        double weight = 1.0 - (**aln).similarity();
        //cout<<"weiht "<<weight<<"\n";
        for(int i=0;i<s.getLength();i++){
            for(int j=i+1;j<s.getLength();j++){
                T bonus = 0.0;
                for(int k=0;k<(**seq).getLength();k++){
                    for(int L=k+1;L<(**seq).getLength();L++){
                        double aln_ik = (**aln).coincidence(i,k);
                        double aln_jL = (**aln).coincidence(j,L);
                        double pair_kL = (**prob)[k][L];
                        bonus += (T) (aln_ik * aln_jL * pair_kL);
                        //cout<<i<<" "<<j<<" "<<k<<" "<<L<<" "<<(T) (aln_ik * aln_jL * pair_kL)<<"\n";
                    }
                }
                bonuses[i][j] += bonus * weight;
                //cout<<bonuses[i][j]<<"\n";
                bonuses[j][i] = bonuses[i][j];
            }
        }
    }
    //cout<<"before normalization\n";
    //show();
	normalize();
	apply_gamma();
}

template<typename T>
T extrinsic<T>::bonus(const int i, const int j) const {
    return initialized ? bonuses[i][j] : ONE;
}

template<typename T>
void extrinsic<T>::show() const{
    cerr <<"extrinsic information"<<endl;
    for(size_t i=0;i<bonuses.size();i++){
        cerr<<"\t"<<i;
    }
    cerr<<endl;
    for(size_t i=0;i<bonuses.size();i++){
        cerr<<i;
        for(size_t j=0;j<bonuses.size();j++){
            cerr<<"\t"<<bonuses[i][j];
        }
        cerr<<endl;
    }
    cerr<<endl;
}

template<typename T>
T extrinsic<T>::max_value() const{
	T tmp = bonuses[0][0];
    for(size_t i=0;i<bonuses.size();i++){
	for(size_t j=0;j<bonuses.size();j++){
		if(bonuses[i][j] > tmp)
			tmp = bonuses[i][j];
	}
	}
	return tmp;
}

template<typename T>
void extrinsic<T>::elementwise_divide(const T denom){
    for(size_t i=0;i<bonuses.size();i++){
	for(size_t j=0;j<bonuses.size();j++){
		bonuses[i][j] /= denom;
	}
	}
}

template<typename T> T exponentiate(T base, T exponent);
template<> int exponentiate<int>(int base, int exponent){
	//this should never happen
	throw "up";
	return 0;
}
template<> real_t exponentiate<real_t>(real_t base, real_t exponent){
	return pow(base, exponent);
}

template<typename T>
void extrinsic<T>::elementwise_pow(const T exponent){
    for(size_t i=0;i<bonuses.size();i++){
	for(size_t j=0;j<bonuses.size();j++){
		bonuses[i][j] = exponentiate(bonuses[i][j], exponent);
	}
	}
}

template<typename T>
void extrinsic<T>::normalize(){
	elementwise_divide(max_value());
}

template<typename T>
void extrinsic<T>::apply_gamma(){
	elementwise_pow(gamma);
}


template class extrinsic<int>;
template class extrinsic<real_t>;
