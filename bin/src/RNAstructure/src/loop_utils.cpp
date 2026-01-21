#include "loop_utils.h"
#include <algorithm>
using std::vector;
using std::sort;
using std::min;
using std::max;

//put this in a namespace to avoid name conflicts with probscan stuff
//would make sense to unify these types at some point
namespace loop{
    basepair::basepair(int i, int j) : i(i),j(j){}

    bool basepair::operator==(const basepair& other) const{
        return this->i==other.i && this->j==other.j;
    }

    basepair& basepair::operator=(const basepair& other){
        this->i = other.i;
        this->j = other.j;
        return *this;
    }

    ostream& operator<<(ostream& output, const basepair& p) {
            output << "(" <<  p.i << ", " << p.j <<")";
            return output;
    }

    loop::loop(int i, int j) : outer(basepair(i,j)){}

    hairpin::hairpin(int i, int j) : loop(i,j) {}

    ostream& operator<<(ostream& output, const hairpin& h) {
            output << "Hairpin: " << h.outer;
            return output;
    }

    vector<int> hairpin::nucs() const
    {
        std::vector<int> ret;
        for(int k=outer.i+1;k<outer.j;k++){
            ret.push_back(k);
        }
        return ret;
    }

    double hairpin::getProbability(ProbScan& p) const{
        return p.probability_of_hairpin(outer.i,outer.j);
    }

    internal::internal(int i, int j, int k, int l) : loop(i,j), inner(k,l) {}

    vector<int> internal::nucs() const
    {
        std::vector<int> ret;
        for(int k=outer.i+1;k<inner.i;k++){
            ret.push_back(k);
        }
        for(int k=inner.j+1;k<outer.j;k++){
            ret.push_back(k);
        }
        return ret;
    }

    double internal::getProbability(ProbScan& p) const{
        return p.probability_of_internal_loop(outer.i,outer.j,inner.i,inner.j);
    }

    ostream& operator<<(ostream& output, const internal& i) {
            output << "Internal: " << i.outer<<" "<<i.inner;
            return output;
    }

    multibranch::multibranch(vector<basepair> pairs) :
        loop(pairs[0].i,pairs[0].j), pairs(pairs) {}

    bool paircomp(const basepair a, const basepair b)
    {
        return a.i < b.i;
    }

    vector<int> multibranch::nucs() const
    {
        vector<int> ret;
        vector<basepair> m = pairs;//non-const working copy
        sort(m.begin(),m.end(),paircomp);
        m[0] = basepair(m[0].j,m[0].i);
        m.push_back(m[0]);
        basepair prev = m[0];
        for(std::vector<basepair>::const_iterator it=m.begin()+1;it!=m.end();++it){
            for(int k = prev.j+1;k<it->i;k++){
                ret.push_back(k);
            }
            prev = *it;
        }
        return ret;
    }

    ostream& operator<<(ostream& output, const multibranch& m) {
            output << "Multibranch: ";
            for(vector<basepair>::const_iterator it=m.pairs.begin();
                    it!=m.pairs.end();
                    ++it){
                output << *it<<" ";
            }
            return output;
    }

    double multibranch::getProbability(ProbScan& p) const
    {
        multibranch_loop_t mb;
        for(vector<basepair>::const_iterator it=pairs.begin();it!=pairs.end();++it){
            std::pair<int,int> p = std::make_pair(it->i,it->j);
            mb.branches.push_back(p);
        }
        return p.probability_of_multibranch_loop(mb);
    }

    stem::stem(int i, int j, int k, int l) : loop(i,j), inner(k,l) {}

    vector<int> stem::nucs() const
    {
        std::vector<int> ret;
        for(int k=outer.i;k<=inner.i;k++){
            ret.push_back(k);
        }
        for(int k=inner.j;k<=outer.j;k++){
            ret.push_back(k);
        }
        return ret;
    }

    double stem::getProbability(ProbScan& p) const{
        return p.probability_of_helix(outer.i,outer.j,inner.i-outer.i);
    }

    ostream& operator<<(ostream& output, const stem& s) {
            output << "Stem: " << s.outer<<" "<<s.inner;
            return output;
    }

    bool unpaired_between(const int i, const int j, RNA& r, const int s){
        for(int k=i+1;k<j;k++){
            if(r.GetPair(k,s)!=0){
                return false;
            }
        }
        return true;
    }

    basepair next_pair(const int i, RNA& r,const int s){
        for(int k = i+1;k<r.GetPair(i,s);k++){
            if(r.GetPair(k,s)!=0){
                return basepair(k,r.GetPair(k,s));
            }
        }
        return(basepair(i,r.GetPair(i,s)));
    }

    bool contains(const basepair outer, const basepair inner){
        return outer.i<inner.i && outer.j > inner.j;
    }

    bool forms_iloop(const basepair outer, const basepair inner, RNA& r, const int s){
        if(outer==inner) return false;
        if(!contains(outer,inner)) return false;
        if(outer.i==inner.i-1 && outer.j==inner.j+1) return false;
        return unpaired_between(outer.i,inner.i,r,s)
            && unpaired_between(inner.j,outer.j,r,s);
    }

    vector<hairpin> find_hairpins(RNA& r,int structurenumber){
        const int n = r.GetSequenceLength();
        const int s = structurenumber;
        vector<hairpin> ret;
        for(int i=1;i<=n;i++){
            if((r.GetPair(i,s)>i) && (unpaired_between(i,r.GetPair(i,s),r,s))){
                ret.push_back(hairpin(i,r.GetPair(i,s)));
            }
        }
        return ret;
    }

    vector<internal> find_internals(RNA& r,int structurenumber){
        const int n = r.GetSequenceLength();
        const int s = structurenumber;
        vector<internal> ret;
        for(int i=1;i<=n;i++){
            if(r.GetPair(i,s)>i){
                const basepair outer = basepair(i,r.GetPair(i,s));
                const basepair inner = next_pair(i,r,s);
                if(forms_iloop(outer,inner,r,s)){
                    ret.push_back(internal(outer.i,outer.j,inner.i,inner.j));
                }
            }
        }
        return ret;
    }

    bool closes_multibranch(const basepair p, RNA& r, const int s){
        int k = p.i+1;
        int count = 0;
        while(k<p.j){
            const int j = r.GetPair(k,s);
            if((j!=0)&&(j<p.i)) return false;
            if(j>p.j) return false;
            if(j>k){
                count +=1;
                k = j;
            }
            else{
                k+=1;
            }
            if(count>10000) {
                std::cerr<<"infinite loop detected\n";
                return false;
            }
        }
        return count > 1;
    }

    multibranch mb_closed_by(const basepair p,RNA& r, const int s){
        std::vector<basepair> pairs;
        pairs.push_back(p);
        int k = p.i+1;
        while(k<p.j){
            if(r.GetPair(k,s)>k){
                pairs.push_back(basepair(k,r.GetPair(k,s)));
                k = r.GetPair(k,s);
            }
            else{
                k+=1;
            }
        }
        return multibranch(pairs);
    }

    vector<multibranch> find_multibranch(RNA& r,const int structurenumber){
        const int s = structurenumber;
        vector<multibranch> ret;
        for(int i=1;i<=(int)r.GetSequenceLength();i++){
            if(r.GetPair(i,s)>i){
                const basepair outer = basepair(i,r.GetPair(i,s));
                const basepair inner = next_pair(i,r,s);
                if(contains(outer,inner) && !forms_iloop(outer,inner,r,s)){
                    if(closes_multibranch(outer,r,s)){
                        ret.push_back(mb_closed_by(outer,r,s));
                    }
                }
            }
        }
        return ret;
    }

    //a nucleotide at position i begins a stem if it is paired to a nucleotide
    //3' of it and either it is the first nucleotide in the sequence or the
    //nucleotide at position i-1 is unpaired
    bool begins_stem(const int i, RNA& r, const int s){
        const int j = r.GetPair(i,s);
        if(j==0 || j<i) return false;
        if(i==1) return true;
        const int jp = r.GetPair(i-1,s);
        if(jp!=j+1 && r.GetPair(i+1,s)==j-1) return true;
        return false;
    }

    //a nucleotide at position i ends a stem if it is paired to a nucleotide
    //3' of it and the nucleotide at position j+1 is unpaired
    bool ends_stem(int i,RNA& r,const int s){
        const int j = r.GetPair(i,s);
        if(j==0 || j<i) return false;
        const int jp = r.GetPair(i+1,s);
        if(jp!=j-1) return true;
        return false;
    }

    vector<stem> find_stems(RNA& r,const int structurenumber){
        const int s = structurenumber;
        vector<stem> ret;
        for(int i=1;i<=r.GetSequenceLength();i++){
            if(begins_stem(i,r,s)){
                int k=i+1;
                while(!ends_stem(k,r,s)){
                    k+=1;
                }
                ret.push_back(stem(i,r.GetPair(i,s),k,r.GetPair(k,s)));
            }
        }
        return ret;
    }

}//end namespace loop
