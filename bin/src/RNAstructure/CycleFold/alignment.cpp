#include "alignment.h"
#include "../src/phmm/phmm_aln.h"
#include <cmath>
#include <iostream>
using std::vector;
using std::cout;
using std::cerr;

vector<char> toVector(const sequence& s){
	vector<char> v;
	for(char c : s.toString()){
		v.push_back(c);
	}
	return v;
}

alignment::alignment(const sequence& a, const sequence& b) :
    probs(vector<vector<double> >(a.getLength(), vector<double>(b.getLength(),0.0)))
{
//perform alignment with forward-backward algorithm
    vector<char> nucs_a = toVector(a);
    vector<char> nucs_b = toVector(b);
    t_phmm_aln* phmm_aln = create_phmm_aln(nucs_a, nucs_b);
    t_pp_result* posterior_probs = phmm_aln->compute_posterior_probs();

//save the sequence similarity
    _similarity = posterior_probs->ml_similarity;

//save co-incidence probabilities in probs
    for(int i=0; i<a.getLength(); i++) {
        for(int k=0; k<b.getLength(); k++) {
          //probabilities are stored as LOGS in the "t_pp_result"
          //so the exponential of the array value is the probability
            double aln_prob = exp(posterior_probs->aln_probs[i][k]);
            double ins1_prob = exp(posterior_probs->ins1_probs[i][k]);
            double ins2_prob = exp(posterior_probs->ins2_probs[i][k]);
            probs[i][k] = aln_prob + ins1_prob + ins2_prob;
        }
    }
	phmm_aln->free_pp_result(posterior_probs);
	delete_phmm_aln(phmm_aln);
  // cout<<"alignment between "<<a.getTag()<<" and "<<b.getTag()<<"\n";
  //show();
}

void alignment::show() const{
  cout<<"similarity "<<_similarity<<"\n";
  cout <<"alignment probabilities"<<endl;
  for(size_t i=0;i<probs.size();i++){
    cout<<"\t"<<i;
  }
  cout<<endl;
  for(size_t i=0;i<probs.size();i++){
    cout<<i;
    for(size_t j=0;j<probs[0].size();j++){
      cout<<"\t"<<probs[i][j];
    }
    cout<<endl;
  }
  cout<<endl;
}


double alignment::coincidence(const int i, const int k) const {
    return probs[i][k];
}

double alignment::similarity() const {
    return _similarity;
}
