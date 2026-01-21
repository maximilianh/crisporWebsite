#include "constraints.h"
using std::vector;
using std::pair;
using std::cout;

vector<int> vrange(int min, int max){
	vector<int> ret;
	for(int i=min;i<max;i++)
		ret.push_back(i);
	return ret;
}

void constraints::show(){
    cout<<"constraints:\n";
    for (size_t i=0;i<pairs.size();i++){
        if(!unpaired(i)){
            cout<<i<<"\t"<<pairs[i]<<"\n";
        }
    }
}

void constraints::show_force_unpaired(){
    cout<<"force unpairing constraints:\n";
    for (size_t i=0;i<unpairs.size();i++){
        if(force_unpaired(i)){
            cout<<i<<"\n";
        }
    }
}

void constraints::build_allowedpairs(const int sequencelength) {

	allowedpairs.resize(sequencelength);
	for (int i=0;i<sequencelength;++i) {
		allowedpairs[i].resize(sequencelength);
		for (int j=0;j<sequencelength;++j) {

			//initialize to true, and subtract any conflicts
			allowedpairs[i][j]=true;

		}
	}


	//Now identify pseudoknot conflicts using array of specified pairs:
	for (int i=0;i<sequencelength;++i) {

		//cout << i << " "<< partner(i)<<"\n";

		if (partner(i)>i) {
			//i is forced paired to partner(i)
			for(int p5=0;p5<i;++p5) {

				

				for (int p3=i+1;p3<partner(i);++p3) {

					//cout << "p5 ="<<p5<< " p3="<<p3<<"\n";

					allowedpairs[p5][p3]=false;
					allowedpairs[p3][p5]=false;
				}

			}
			for(int p5=i+1;p5<partner(i);++p5) {

				

				for (int p3=partner(i)+1;p3<sequencelength;++p3) {

					//cout << "p5 ="<<p5<< " p3="<<p3<<"\n";

					allowedpairs[p5][p3]=false;
					allowedpairs[p3][p5]=false;
				}

			}

		}
	}



}

bool constraints::IsPairAllowed(const int i, const int j) const{

	return (!constrained) || allowedpairs[i][j];

}

vector<int> unpaired_indices(const int i, const int j, const NCM_type ncm){
    vector<int> ret;
	int start_five = i<j? i : j-NCM::fivep_length(ncm)+1;
    for(int k=start_five+1; k<start_five+NCM::fivep_length(ncm)-1; k++){
        ret.push_back(k);
    }
	int start_three = i<j? j : i+NCM::threep_length(ncm)-1;
    for(int k=start_three-1; k>start_three-NCM::threep_length(ncm)+1; k--){
        ret.push_back(k);
    }
    return ret;
}

vector<int> unpaired_indices(const int i, const int j, const int k, const int l){
    vector<int> ret;
	bool interior = i<j;
	int start_five = interior? i:l;
	int end_five = interior? k:j;
    for(int kp=start_five+1; kp<end_five; kp++){
        ret.push_back(kp);
    }
	int start_three = interior? l:i;
	int end_three = interior? j:k;
    for(int kp=start_three+1; kp<end_three; kp++){
        ret.push_back(kp);
    }
    return ret;
}

vector<int> paired_indices(const int i, const int j, const NCM_type ncm){
    vector<int> ret = {i, j};
    if(ncm[0]=='2'){
		if(i<j){
			ret.push_back(i+NCM::fivep_length(ncm)-1);
			ret.push_back(j-NCM::threep_length(ncm)+1);
		}
		else {
			ret.push_back(i+NCM::threep_length(ncm)-1);
			ret.push_back(j-NCM::fivep_length(ncm)+1);
		}
    }
    return ret;
}

bool conflicts(const int i,const int j, const int k, const int l, const constraints& cst){
    if(!cst.unpaired(i) || !cst.unpaired(j)){ //if nucs at i or j  are paired to something..
        if(cst.partner(i) != j){ //and that something is not j..
            return true; //there is a conflict
        }
    }
	if (!cst.IsPairAllowed(i,j)) return true;
    for(int kp : unpaired_indices(i,j,k,l)){ //if ncm at i,j implies k is unpaired..
        if(cst.valid(kp) && !cst.unpaired(kp)){ //but k is actually paired to something..
            return true; //there is a conflict
        }
	}
	return false;
}

bool conflicts(const int i, const int j, const NCM_type ncm, const constraints& cst){
    if(!cst.unpaired(i) || !cst.unpaired(j)){ //if nucs at i or j  are paired to something..
        if(cst.partner(i) != j){ //and that something is not j..
            return true; //there is a conflict
        }
    }
	if (!cst.IsPairAllowed(i,j)) return true;
    for(int k : unpaired_indices(i,j,ncm)){ //if ncm at i,j implies k is unpaired..
        if(cst.valid(k) && !cst.unpaired(k)){ //but k is actually paired to something..
            return true; //there is a conflict
        }
    }
    for(int k : paired_indices(i,j,ncm)){
        if(cst.valid(k) && cst.force_unpaired(k)){
            return true;
        }
    }
    return false; //anything that is not forbidden is allowed
}


constraints::constraints(): constrained(false) {}

constraints::constraints(vector<int> p):
	constrained(true),
	pairs(p),
	unpairs(vector<bool>(p.size(),false)) {}

constraints::constraints(vector<bool> u):
	constrained(true),
	pairs(vrange(0,u.size())),
	unpairs(u) {}

bool constraints::valid(const int i) const{
	return i>=0 && i<(int)pairs.size();
}

bool constraints::force_unpaired(const int i) const {
    return constrained && unpairs[i];
}

bool allowed_pair(const int i, const int j, const constraints& cst){
    if(cst.force_unpaired(i) || cst.force_unpaired(j)){
        return false;
    }
    return cst.partner(i) == j || (cst.unpaired(i) && cst.unpaired(j));
}

bool allowed_unpair(const int i, const constraints& cst){
    return cst.unpaired(i);
}

bool constraints::unpaired(const int i) const {
    if(constrained){
        return pairs[i] == i;
    }
    return true;
}

int constraints::partner(const int i) const {
    if(constrained){
        return pairs[i];
    }
    return i;
}
