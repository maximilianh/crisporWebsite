#include "turbo_calculation.h"
#include "alignment.h"
#include "mainloop.h"
#include "extrinsic.h"
#include "constraints.h"
#include <map>

using std::map;
using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;


map<int,vector<alignment>> setup_alignments(const vector<sequence>& seqs){
    map<int,vector<alignment>> alignments;
    for (size_t i=0;i<seqs.size();i++){
        vector<alignment> aln;
        for (size_t j=0;j<seqs.size();j++){
            if(i==j) continue;
            //cout<<"aignment between "<<i<<" and "<<j<<"\n";
            //alignment(seqs[i],seqs[j]).show();
            aln.push_back(alignment(seqs[i],seqs[j]));
        }
        alignments[i] = aln;
    }
    return alignments;
}

//This builds out the table_t pair probabilities by allocating space and by running the CycleFold partition function.
vector<table_t> calculate_initial_pair_probs(const vector<sequence>& seqs,
                                             const parameters<real_t>& p){
    vector<table_t> initial_probs;
    for(const sequence& seq : seqs){
        initial_probs.push_back(calculate_pairing_probabilities(seq,p,constraints()));
    }
    return initial_probs;
}

template <typename T>
vector<const T*> other(const T& x, const vector<T>& xs){
    vector<const T*> others;
    for(const T& other: xs){
        if(&other != &x){
            others.push_back(&other);
        }
    }
    return others;
}

//return a new vector with pointers to everything in xs
//except the item at index x
template <typename T>
vector<const T*> other(const size_t x, const vector<T>& xs){
  vector<const T*> others;
  for(size_t i=0; i<xs.size(); i++){
    if(i!=x){
      others.push_back(&(xs[i]));
    }
  }
  return others;
}

template <typename T>
vector<const T*> pointerfy(const vector<T>& xs){
    vector<const T*> ret;
    for(const T& x : xs){
        ret.push_back(&x);
    }
    return ret;
}

//Return a vector (dimension is sequence number) of table_t pair probabilities (each is 2D array).
vector<table_t> turbo_calculation(const vector<sequence>& seqs,
                                  const parameters<real_t>& p,
                                  const int n_iterations,
                                  const real_t gamma,
                                  const vector<constraints>& constr)
{
	//Make a vector of pairwise alignment probabilities.  This is the original TurboFold style, which does not refine
	//the alignment probabilities.
    const map<int,vector<alignment>> alignments = setup_alignments(seqs);
	
	//This runs the CycleFold partition function for each sequence.
    vector<table_t> probs = calculate_initial_pair_probs(seqs,p);

	//Refine pair probabilities for each iteration.
    for(int iter=0; iter<n_iterations; iter++){
        vector<table_t> new_probs(seqs.size());
        for(int i=0; i<(int)seqs.size(); i++){

			//use pointers to avoid having to copy the sequences, probs
			//and alignments at every iteration
            const vector<const sequence*> other_seqs = other(i, seqs);
            const vector<const table_t*> other_probs = other(i, probs);
            const vector<const alignment*> other_alns = pointerfy(alignments.at(i));
            const extrinsic<real_t> ext = extrinsic<real_t>(seqs[i], other_seqs, other_probs, other_alns,gamma);
            //ext.show();
            constraints c = constr.empty()? constraints() : constr[i];
            table_t new_prob = calculate_pairing_probabilities(seqs[i], p, c, ext);
            new_probs[i] = new_prob;
            //cerr<<"probs\n";
            //show_pair_probs_arrayview(new_prob);
        }

		//update the probabilities stored in probs with the updated values.
        probs = new_probs;
    }
    //for n_iterations
    //  for n_seqs
    //      calculate extrinsic info with old probs
    //      calculate new probs
    //  update old probs with new probs

    return probs;
}
