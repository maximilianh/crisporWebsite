#include <math.h>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "domain_class.h"


using namespace std;

domain_calc::domain_calc(bool circularize, double outside_weight, int leng) {
	length=leng;
	//min_domain_length = 20;

	//Initialize DynProgArrayU objects
	prob_matrix = new DynProgArrayU<double>(length, 0);
	row_sum = new DynProgArrayU<double>(length, 0);
	col_sum = new DynProgArrayU<double>(length, 0);
	inner = new DynProgArrayU<double>(length, 0);
	total = new DynProgArrayU<double>(length, 0);
	scores = new DynProgArrayU<double>(length, -1);
	
	if (circ)
		scores2 = new DynProgArrayU<double>(length, -1);
/*
	prob_matrix = new DynProgArrayU<double>(length, 0);
	row_sum = new DynProgArrayU<double>(length, 0);
	col_sum = new DynProgArrayU<double>(length, 0);
	inner = new DynProgArrayU<double>(length, 0);
	total = new DynProgArrayU<double>(length, 0);
	scores = new DynProgArrayU<double>(length, -1);
*/
	circ = circularize;

	w2 = outside_weight;

}

domain_calc::~domain_calc() {
//	delete circ;

//	delete w2;

//	delete base_pairs;

	//The sequence of the RNA to be segmented
//	delete seq;

	//The length of the RNA, in nucleotides
//	delete length;

	//A number of dynamic programming arrays used by the domain findng algorithm
	delete prob_matrix;
	delete row_sum;
	delete col_sum;
	delete inner;
	delete total;
	delete scores;
	
	//constraints contains a number of i,j pairs, where (i,j) defines a base pair closing a folding domain
//	delete constraints;
}


//pfunction performs the partition function calculation and saves the results to disk.
void domain_calc::fill_domain_matrices(){
/* 	row_sum = new domainArray(length);
	col_sum = new domainArray(length);
	inner = new domainArray(length);
	total = new domainArray(length);
	prob_matrix = new domainArray(length);
	scores = new domainArray(length,-1);
 */
	int i, j, k, n;

	//Add base pair probabilities to the probability matrix
	for(i=0; i < base_pairs.size(); i++){
		prob_matrix->f(base_pairs[i].i,base_pairs[i].j)=base_pairs[i].prob;
	}

	if (circ){
		int max_n = (int) length/2+0.5;
	}
	else{
		int max_n = length;
	}

	for(n=0; n < length; n++){
		for(i=0; i < length-n-1; i++){
			j = i+n+1;
			row_sum->f(i,j)=row_sum->f(i,j-1)+prob_matrix->f(i,j);
			col_sum->f(i,j)=col_sum->f(i+1,j)+prob_matrix->f(i,j);
		}
	}

	for (i=0; i<length; i++){
        total->f(i,i) = row_sum->f(i,length-1)+col_sum->f(0,i);
	}

	for(n=0; n < length; n++){
		for(i=0;i < length-n-1;i++){
			j = i+n+1;

			inner->f(i,j)=inner->f(i+1,j-1)+row_sum->f(i,j)+col_sum->f(i,j)-prob_matrix->f(i,j);
            total->f(i,j)=total->f(i+1,j-1)+row_sum->f(i,length-1)-row_sum->f(i,j)+col_sum->f(0,i)+row_sum->f(j,length-1)+col_sum->f(0,j)-col_sum->f(i+1,j);
		}
	}

}

//writepfsave writes a save file with partition function data.
void domain_calc::writeconstraints(std::string filename){
	ofstream ofile;
	ofile.open(filename.c_str());

	ofile << "# 5'\t3'\tScore" << endl;

	for(int i=0; i<constraints5.size(); i++){
		ofile << constraints5[i]+1 << "\t" << constraints3[i]+1 << "\t" << score(constraints5[i], constraints3[i]) << endl;
	}//*/
}

//readpfsave reads a save file with partition function data.
void domain_calc::read_pp_file(std::string filename){
	pair_data temp_pair_info;
	
	int i, j;
	double prob;
	
	ifstream pp_file;
	pp_file.open(filename.c_str(), ios::in);

	//ifstream pp_file(filename);

	string line;

	//skip first two lines
	for(i=0; i<2; i++){
		getline(pp_file,line);
	}

	while(pp_file >> i >> j >> prob){
		if (j > i) {
			temp_pair_info.i = i-1;
			temp_pair_info.j = j-1;
			temp_pair_info.prob = prob;

			base_pairs.push_back(temp_pair_info);
		}
	}

	//return 0
}

void domain_calc::read_pfs_file(RNA *strand){
	pair_data temp_pair_info;

	seq = strand->GetSequence();
	length = strand->GetSequenceLength();

	// Add all possible pairs
	for( int i = 1; i <= length; i++ ) {
		for( int j = i+1; j <= length; j++ ) {
			double value = strand->GetPairProbability( i, j );
			if( ( value != 0.0 ) && ( value != -0.0 ) ){
				temp_pair_info.i = i-1;
				temp_pair_info.j = j-1;
				temp_pair_info.prob = value;

				base_pairs.push_back(temp_pair_info);
			}
		}
	}
}

void domain_calc::calc_domain_scores(){
	int j;

	for(int n=0; n < length; n++){
		for(int i=0; i < length-n; i++){
			j = i+n;

			// For circular domains:
			//		inner = total(0,length-1)-total(i+1,j-1)
			//		outer = total(i+1,j-1)-inner(i+1,j-1)
			//		total = total(0,length - inner(i+1,j-1)
			if (circ)
				scores2->f(i,j) = (total->f(0,length-1)-(1+w2)*total->f(i+1,j-1)+w2*inner->f(i+1,j-1))/(total->f(0,length-1)-inner->f(i+1,j-1));
//				scores2->f(i,j) = (total->f(0,length-1)-(2+w2)*total->f(i+1,j-1)+(1+w2)*inner->f(i+1,j-1))/(total->f(0,length-1)-total->f(i+1,j-1));
            
			//score = (inner - w2*outer)/total = (inner-w2*(total-inner))/total
			scores->f(i,j) = ((1+w2)*inner->f(i,j)-w2*total->f(i,j))/(total->f(i,j));
		}
	}
}

void domain_calc::complete_domain_assignment(int min_domain_length){
	int max_phase;
	int d;

	//cout << "Determine max_phase" << endl;
	// Determine the maximum number of phase
	if (circ)
		max_phase = length/2;
	else
		max_phase = 1;
	
	// Declare the optimization array
	//cout << "Declare DP arrays" << endl;
	double** V_array = new double*[max_phase];//[max_phase][length];
	int** Domain_count = new int*[max_phase];//[max_phase][length];	

	for (int i=0; i<max_phase; i++){
		V_array[i] = new double[length+1];
		//Domain_count[i] = new int[length+1];
	}


	//cout << "Initialize Optimization Array" << endl;
	// Initialize the optimization array
	for(d=0; d<max_phase; d++){
		V_array[d][0]=0;			
		//Domain_count[d][0]=0;
		for (int i=1; i<length+1; i++){
			V_array[d][i]=-10000;			
			//Domain_count[d][i]=0;
		}
	}
	//cout << "Fill Array" << endl;
	// Fill the optimization array
	for(d=0; d<max_phase; d++){
		//cout << "Phase" << d << endl;
		for (int i=min_domain_length+1; i<=length; i++){
			for (int m=1; m<=i-min_domain_length; m++){
				
//				if (V_array[d][i] < pow(sqrt(V_array[d][m]*Domain_count[d][m])+score(m+d,i+d),2)/(Domain_count[d][m]+1)){
//					V_array[d][i] = pow(sqrt(V_array[d][m]*Domain_count[d][m])+score(m+d,i+d),2)/(Domain_count[d][m]+1);
//					Domain_count[d][i] = Domain_count[d][m]+1;
//				}

//				if (V_array[d][i] < V_array[d][m]+pow(score(m+d,i+d),3)){
//					V_array[d][i] = V_array[d][m]+pow(score(m+d,i+d),3);
//					Domain_count[d][i] = Domain_count[d][m]+1;
//				}
				if (V_array[d][i] < V_array[d][m-1]+score(m+d-1,i+d-1)){
					V_array[d][i] = V_array[d][m-1]+score(m+d-1,i+d-1);
					//if (d == 21)
					//	cout << V_array[d][m] << "\t" << score(m+d+1,i+d) << V_array[d][i] << endl;
				}
			}
		}
	}

	domain_phase = 0;
	//cout << "Find Optimum phase" << endl;
	if (circ){
		for(int i=0;i<max_phase; i++){
			//cout << i << ":   " << V_array[i][length-1] << endl;
			if (V_array[i][length]>V_array[domain_phase][length])
				domain_phase = i;
		}
	}
	
	d = domain_phase;

	cout << "Max Phase:    " << d << endl;
	cout << "Max Score:    " << V_array[d][length] << endl;
	cout << "Length:       " << length << endl;
	//cout << "Domain Count: " << Domain_count[d][length] << endl;

	// Start traceback
	// Find the optimum phase
	//cout << "Start Traceback" << endl;	
	
	int i  = length;
	bool loop_flag = false;

	int domain5, domain3;

	while (i != 0 && !loop_flag){
        loop_flag = true;
		for (int m = 1; m<=i-min_domain_length; m++){
            if (V_array[d][i] == V_array[d][m-1]+score(m+d-1,i+d-1)){

//            if (V_array[d][i] == pow(sqrt(V_array[d][m]*Domain_count[d][m])+score(m+d,i+d),2)/(Domain_count[d][m]+1)){
//            if (V_array[d][i] == V_array[d][m]+pow(score(m+d,i+d),3)){
				domain5 = m+d;
				domain3 = i+d;

				if (domain5 > length)
					domain5 -= length;
				if (domain3 > length)
					domain3 -= length;

				//if (score(m+d,i+d)>0){
					if (domain3 > domain5){
						constraints5.push_back(domain5-1);
						constraints3.push_back(domain3-1);
					}
				//}
				cout << domain5 << "\t";
				cout << domain3 << "\t";
				
				if (domain5 < domain3)
					cout << score(m+d-1,i+d-1) << "\t" << inner->f(m+d-1,i+d-1) <<"\t" << total->f(m+d-1,i+d-1) << endl;
				else
					cout << score(m+d-1,i+d-1) << "\t" << total->f(0,length-1)-total->f(domain3+1,domain5-1) <<"\t" << (total->f(0,length-1)-inner->f(domain3+1,domain5-1)) << endl;
			
                i = m-1;
				loop_flag = false;
                break;
			}
		}		
		if (loop_flag){
			cout << "Infinite loop detected" << endl;
			//for(int m = 0;m<=i;m++){
			//	cout << m << "\t" << V_array[d][m] << endl;
			//}
			constraints5.push_back(d);
			constraints3.push_back(i+d);
			cout << d+1 << "\t";
			cout << i+d+1 << "\t";
			cout << pow(score(d,i+d),1) << endl;
		}
	}
}

double domain_calc::score(int i, int j){
	int i_eff, j_eff;

	if (i> length-1)
        i_eff = i-length;
    else
        i_eff = i;
        
    if (j > length-1)
        j_eff = j-length;
    else
        j_eff = j;
        
    if (j_eff < i_eff){
		if (circ)
			return scores2->f(j_eff, i_eff);
		else{
			cout << "Bad access" << endl;
			return -10000;
		}
    }
	else
		return scores->f(i_eff, j_eff);
}

void domain_calc::copy_constraints(structure *ct){
	for (int i=0; i<constraints5.size(); i++){
		ct->AddDomain(constraints5[i]+1, constraints3[i]+1);
	}
}
