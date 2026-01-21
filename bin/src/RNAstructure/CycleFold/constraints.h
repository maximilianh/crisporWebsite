#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H
#include "NCM_parameters.h"
#include <vector>

//an object that stores secondary structrue constraint information
class constraints{
    public:
        constraints();
        constraints(std::vector<int> p);
        constraints(std::vector<bool> u);
        constraints(std::vector<std::pair<int,int>> p);
		bool valid(const int i) const;
        bool unpaired(const int i) const;
        bool force_unpaired(const int i) const;
        int partner(const int i) const;
        void show();
        void show_force_unpaired();
		void build_allowedpairs(const int sequencelength);//build the allowedpairs array, requires the length of the sequence
		bool IsPairAllowed(const int i, const int j) const;//Return true if i-j is allowed, given the set of pair constraints 
    private:
        bool constrained;
        std::vector<int> pairs;
        std::vector<bool> unpairs;
		std::vector< std::vector<bool> > allowedpairs;//this new array will track pseudoknot conflicts with forced pairs
};

//given an NCM at a position, returns the list of unpaired nucleotides in that NCM
std::vector<int> unpaired_indices(const int i, const int j, const NCM_type ncm);
//given an NCM at a position, returns the list of paired nucleotides in that NCM
std::vector<int> paired_indices(const int i, const int j, const NCM_type ncm);

bool allowed_pair(const int i, const int j, const constraints& cst);
bool allowed_unpair(const int i, const constraints& cst);

//given an NCM at a position, returns whether that NCM would require any forced pairs to be unpaired
//i.e. this structure is impossible
bool conflicts(const int i, const int j, const NCM_type ncm, const constraints& cst);

#endif
