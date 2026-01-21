#include "maxexpect.h"
#include <vector>
#include <algorithm>
#include <stack>
#include <cmath>
#include <limits>
const bool VERBOSE = false;
using std::vector;
using std::max;
using std::pair;
using std::stack;
using std::cout;
using std::endl;

void show(const table_t& t){
    for(size_t i=0;i<t.size();i++){
        for(size_t j=i;j<t.size();j++){
            cout<<"i: "<<i<<" j: "<<j<<" value: "<<t[i][j]<<endl;
        }
    }
}

vector<double> unpairing_probabilities(const table_t& pair_probs){
    vector<double> q(pair_probs.size());
    for(size_t i=0;i<pair_probs.size();i++){
        double p = 1.0;
        for(size_t j=0;j<pair_probs.size();j++){
            if(i==j) continue;
            p -= pair_probs[i][j];
        }
        q[i] = p;
    }
    return q;
}

table_t maxexpect_fill(const table_t& p, const vector<double>& q, double gamma){
    const int n = p.size();
    vector<vector<double>> M(n,vector<double>(n,0.0));
    for(int i=n-1;i>=0;i--){
        for(int j=i;j<n;j++){
            double tmp = 0.0;
            if(i==j){
                tmp = max(tmp,q.at(i));
            }
            if(i<j && i<n-1){
                tmp = max(tmp,q.at(i)+M.at(i+1).at(j));
            }
            if(i<j && j>0){
                tmp = max(tmp,q.at(j)+M.at(i).at(j-1));
            }
            if(i<j && j>0 && i<n-1){
                tmp = max(tmp,gamma*2*p.at(i).at(j)+M.at(i+1).at(j-1));
            }
            for(int k=i;k<j;k++){
                tmp = max(tmp,M.at(i).at(k)+M.at(k+1).at(j));
            }
            M[i][j] = tmp;
        }
    }
    if(VERBOSE) show(M);
    return M;
}

bool CLOSE(const double a, const double b){
    return fabs(a-b) < 100*std::numeric_limits<double>::min();
}

int branch(const int i, const int j, const table_t& M){
    const double mij = M.at(i).at(j);
    for(int k=i;k<j;k++){
        if(CLOSE(mij, M.at(i).at(k)+M.at(k+1).at(j))){
            return k;
        }
    }
    return -1;
}

vector<pair<int,int>> maxexpect_traceback(const table_t& M, const table_t& p, const vector<double>& q, const double gamma){
    const int n = p.size();
    vector<pair<int,int>> pairs;
    stack<pair<int,int>> stk;
    stk.push({0,n-1});
    while(!stk.empty()){
        pair<int,int> pr = stk.top();
        const int i = pr.first;
        const int j = pr.second;
        int k = 0; //position of a branch if this is a bifurcation
        const double mij = M.at(i).at(j);
        stk.pop();

        //check the possible ways we could have reached this array value
        //hairpin
        if(j-i < 2){
            if (VERBOSE) cout<<"hairpin at "<<i<<" "<<j<<endl;
            continue;
        }
        //base pair
        else if (CLOSE(mij, gamma*2*p.at(i).at(j)+M.at(i+1).at(j-1))){
            if (VERBOSE) cout << "found a pair at "<<i<<" "<<j<<endl;
            pairs.push_back({i,j});
            stk.push({i+1,j-1});
        }
        //5' neighbor
        else if (CLOSE(mij, q.at(i)+M.at(i+1).at(j))){
            if (VERBOSE) cout<<"5' neighbor at "<<i<<" "<<j<<endl;
            stk.push({i+1,j});
        }
        //3' neighbor
        else if (CLOSE(mij, q.at(j)+M.at(i).at(j-1))){
            if (VERBOSE) cout<<"3' neighbor at "<<i<<" "<<j<<endl;
            stk.push({i,j-1});
        }
        //bifurcation
        else if ( (k = branch(i,j,M)) != -1 ){
            if (VERBOSE) cout<<"branch "<<k<<" at "<<i<<" "<<j<<endl;
            stk.push({i,k});
            stk.push({k+1,j});
        }
        else {
            throw "error in maxexpect traceback\n";
        }
    }
    return pairs;
}

vector<pair<int,int>> maxexpect(const table_t& p, double gamma /* =1.0 */){
    const vector<double> q = unpairing_probabilities(p);
    const table_t M = maxexpect_fill(p,q,gamma);
    const vector<pair<int,int>> structure = maxexpect_traceback(M,p,q,gamma);
    return structure;
}
