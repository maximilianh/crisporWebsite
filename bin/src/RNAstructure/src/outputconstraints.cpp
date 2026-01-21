
#include "outputconstraints.h"
#include <stdio.h>
#include <fstream>
#include <cstring>

using namespace std;





void outputconstraints(const char *filename, structure *ct) {
	ofstream out;
	int i,j,ip;

	out.open(filename);
		
	//output ds
	out << "DS:\n";
	for (i=0;i<ct->GetNumberofDoubles();i++) {
		out << ct->GetDouble(i);
		out << "\n";
	}
	out << "-1\n";

	//output ss


	out << "SS:\n";
	for (i=0;i<ct->GetNumberofSingles();i++) {
		out << ct->GetSingle(i);
		out << "\n";
	} 
	
	out << "-1\n";
	

	//output mod
	
	out << "Mod:\n";
	for (i=0;i<ct->GetNumberofModified();i++) {
		out << ct->GetModified(i);
		out << "\n";
	} 

	out << "-1\n";

	//output pairs
	out << "Pairs:\n";
	for (i=0;i<ct->GetNumberofPairs();i++) {
		out << ct->GetPair5(i) << " "<<ct->GetPair3(i);
		out << "\n";
	} 

	out << "-1 -1\n";


	//output fmn

	out << "FMN:\n";
	for (i=0;i<ct->GetNumberofGU();i++) {
		out << ct->GetGUpair(i);
		out << "\n";
	} 

	out << "-1\n";

	//output forbidden pairs
	out << "Forbids:\n";
	for (i=0;i<ct->GetNumberofForbiddenPairs();i++) {
		out << ct->GetForbiddenPair5(i)<< " "<<ct->GetForbiddenPair3(i);
		out << "\n";
	} 

	out << "-1 -1\n";

	//output NMR constraints -- if they exist
	if (ct->min_g_or_u>0||ct->min_gu>0||ct->nneighbors>0||ct->nregion>0) {
		out << "Minimum G or U pairs:\n";
		out << ct->min_g_or_u << "\n";
		out << "Minimum GU pairs:\n";
		out << ct->min_gu << "\n";
		out << "Neighbors:\n";
		for (i=0;i<ct->nneighbors;i++) {
			for (j=0;ct->neighbors[i][j]>0;j++) {
				out << ct->neighbors[i][j] << " ";
			}
			out << "-1\n";

		}
		out << "-1\n";

		//now output the regional nmr constraint data
		out << "Number of NMR Constraint Regions: \n";
		out << ct->nregion<<"\n";
		for (ip=0;ip<ct->nregion;ip++) {
			out << "Start:\n";
			out << ct->start[ip]<<"\n";
			out << "Stop:\n";
			out << ct->stop[ip]<<"\n";
			out << "     Minimum G or U pairs:\n";
			out << ct->rmin_g_or_u[ip] << "\n";
			out << "     Minimum GU pairs:\n";
			out << ct->rmin_gu[ip] << "\n";
			out << "     Neighbors:\n";
			for (i=0;i<ct->rnneighbors[ip];i++) {
				for (j=0;ct->rneighbors[ip][i][j]>0;j++) {
					out << ct->rneighbors[ip][i][j] << " ";
				}
				out << "-1\n";

			}
			out << "-1";

		}

	}

	//output the microarray constraints
	out << "Microarray Constraints:\n";
	out << ct->nmicroarray <<"\n";
	for (ip=0;ip<ct->nmicroarray;ip++) {
		out << ct->microstart[ip] << " " << ct->microstop[ip] << " " << ct->microunpair[ip] << "\n";

	}


	out.close();

}

bool readconstraints(const char *filename, structure *ct) {

	ifstream in;
	int i,ip,constraint,constraint2;

	char temp[40];


	in.open(filename);

	//input ds
	in >> temp;
	

	//Read the double stranded constraints:
	in>> constraint;
	while (constraint!=-1) {
		
		ct->AddDouble(constraint);
		
		in>> constraint;
	}
	

	//input ss
	in >> temp;


	//Read the single-stranded constraints:
	in>> constraint;
	while (constraint!=-1) {
		
		ct->AddSingle(constraint);

		in>> constraint;
	}
	
	

	//input mod
	in >> temp;

	//Read the nucleotides accessible to chemical modification:
	
	in>> constraint;
	while (constraint!=-1) {

		ct->AddModified(constraint);
		
		in>> constraint;
	}
	

	//input pairs
	in >> temp;

	//Read the pairs
	
	in>> constraint;
	in>> constraint2;
	while (constraint!=-1) {
		
		ct->AddPair(constraint,constraint2);
		
		in>> constraint;
		in>> constraint2;

		
	}
	


	//input fmn
	in >> temp;

	//Read the set of Us in GU pairs:
	
	in>> constraint;
	while (constraint!=-1) {
		
		ct->AddGUPair(constraint);
		
		in>> constraint;
	}
	

	//input forbidden pairs
	in >> temp;
	
	//Read the forbidden pairs
	in>> constraint;
	in>> constraint2;
	while (constraint!=-1) {
		
		ct->AddForbiddenPair(constraint,constraint2);
		
		in>> constraint;
		in>> constraint2;
		
	}
	
	//in >> temp;
	//input NMR constraints--if they were written
	in.getline(temp,39);
	in.getline(temp,39);
	if (!in.eof()) {
		if (!strcmp(temp,"Minimum G or U pairs:")) {

			in>> ct->min_g_or_u;
			in>>temp;
		
			in.getline(temp,19);
			in>> ct->min_gu;
			in >>temp;

		
			i = 0;
			ct->nneighbors = 0;
			in>> ct->neighbors[ct->nneighbors][i];
			while (ct->neighbors[ct->nneighbors][i]!=-1) {
				i++;
				in>> ct->neighbors[ct->nneighbors][i];
				while (ct->neighbors[ct->nneighbors][i]!=-1) {
					i++;
					in>> ct->neighbors[ct->nneighbors][i];
				}
				ct->neighbors[ct->nneighbors][i]=0;
				ct->nneighbors++;
				in>> ct->neighbors[ct->nneighbors][i];

			}
			//now input the regional nmr constraint data
			in >>temp;
			in.getline(temp,39);
			in >> ct->nregion;
			for (ip=0;ip<ct->nregion;ip++) {
				in >>temp;
				in >> ct->start[ip];
				in >>temp;
				in >> ct->stop[ip];
				in >>temp;
				in.getline(temp,19);
				in >> ct->rmin_g_or_u[ip];
				in >>temp;
				in.getline(temp,19);
				in >> ct->rmin_gu[ip];
				ct->rnneighbors[ip] = 0;
				in >>temp;
				in>> ct->rneighbors[ip][ct->rnneighbors[ip]][i];
				while (ct->rneighbors[ip][ct->rnneighbors[ip]][i]!=-1) {
					i++;
					in>> ct->rneighbors[ip][ct->rnneighbors[ip]][i];
					while (ct->rneighbors[ip][ct->rnneighbors[ip]][i]!=-1) {
						i++;
						in>> ct->rneighbors[ip][ct->rnneighbors[ip]][i];
					}
					ct->rneighbors[ip][ct->rnneighbors[ip]][i]=0;
					ct->rnneighbors[ip]++;
					in>> ct->rneighbors[ip][ct->rnneighbors[ip]][i];

				}

			
			}
			in.getline(temp,39);
		}
		else ct->min_g_or_u = 0;

	}
	else ct->min_g_or_u = 0;
	if (!in.eof()) {
		if (!strcmp(temp,"Microarray Constraints:")) {
			in >> ct->nmicroarray;
			for (ip=0;ip<ct->nmicroarray;ip++) {
				in >> ct->microstart[ip];
				in >> ct->microstop[ip];
				in >> ct->microunpair[ip];

			}

		}
	}
	else ct->nmicroarray=0;
	


	in.close();
	return true;  //this can be changed later to return error messages
}


