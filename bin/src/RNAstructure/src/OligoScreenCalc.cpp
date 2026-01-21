
#include "OligoScreenCalc.h"
#include "algorithm.h"
#include <string>
#include <cstring>
#include <fstream>

#ifdef SMP
#include <omp.h>
#endif

using namespace std;

//Perform the backend calculation for OligoScreen

//Return an int that indicates an error
//0 = no error

//Pre-define the suboptimal structure prediction parameters.
#define CNTR6 100 //tracebacks
#define CNTR8 20 //percent sort
#define CNTR9 0 //window size


int bimolecular(structure *ct,datatable *data) {
	//Calculate the bimolecular folding free energy
	structure *ct2;
	int i,energy;

	ct2 = new structure();

	//Set the location of the thermodynamic data tables
	ct2->SetThermodynamicDataTable(data);
	

	

	//first create a structure with the sequence needed for bimolecular folding
	ct2->allocate(2*ct->GetSequenceLength()+3);
	ct2->intermolecular = true;
	for (i=1;i<=ct->GetSequenceLength();i++) {
		ct2->numseq[i] = ct->numseq[i];
		ct2->numseq[i+ct->GetSequenceLength()+3] = ct->numseq[i];


	}

	//Set the linker, using the first (and probably only) linker entry

	ct2->numseq[ct->GetSequenceLength()+1] = data->basetonum(data->linker[0]);
	ct2->numseq[ct->GetSequenceLength()+2] = data->basetonum(data->linker[0]);
	ct2->numseq[ct->GetSequenceLength()+3] = data->basetonum(data->linker[0]);

	ct2->inter[0] = ct->GetSequenceLength()+1;
	ct2->inter[1] = ct->GetSequenceLength()+2;
	ct2->inter[2] = ct->GetSequenceLength()+3;

	//ct2->GetSequenceLength() = 2*ct->GetSequenceLength()+3;
	ct2->SetSequenceLabel(ct->GetSequenceLabel());
	//strcpy(ct2->ctlabel[1],ct->ctlabel[1]);
		
	//now fold and perform efn2
	dynamic(ct2,data,CNTR6,CNTR8,CNTR9,NULL,true);
	
	
	energy = ct2->GetEnergy(1);
	delete ct2;
	return energy;

}

//Perform the OligoScreen calculation.
int OligoScreenCalc(const char *infilename, const char *outfilename, datatable *data, rddata *hybriddata) {
	ifstream in;
	ofstream out;
//<<<<<<< OligoScreenCalc.cpp
//	char temp[500];
//	int i,j,energy;
//	structure *ct;


//=======
#ifdef SMP
	int threads = omp_get_num_procs();
#else
	int threads = 1;
#endif
	string *oligo;
	oligo = new string[threads];
//>>>>>>> 1.7
	in.open(infilename);
	out.open(outfilename);
	
	out << "Sequence\tDGbimolecular\tDGunimolecular\tDGduplex\tDG2BPat5'\tDG2BPat3'\n";


	while (!in.eof()) {
		for (int i=0;i<threads;i++){
			if (!in.eof()) {
				getline(in, oligo[i]);
			}
			else 
				oligo[i].erase();
		}

//<<<<<<< OligoScreenCalc.cpp
///		out << temp << "\t";
		//convert temp to a ct file
//		i = strlen(temp);
//		ct = new structure();

//		ct->allocate(i);
		//ct->numofbases=i;

//		ct->SetSequenceLabel("");
		//strcpy(ct->ctlabel[1],"");
//		for (j=0;j<i;j++) {

//			if (temp[j]=='A'||temp[j]=='a') ct->numseq[j+1] = 1;
//			else if (temp[j]=='C'||temp[j]=='c') ct->numseq[j+1] = 2;
//			else if (temp[j]=='G'||temp[j]=='g') ct->numseq[j+1] = 3;
//			else  ct->numseq[j+1] = 4;

//		}
//=======
#ifdef SMP
#pragma omp parallel for ordered schedule(runtime)
#endif
		for (int i=0;i<threads;i++){	
			if (!oligo[i].empty()){
//>>>>>>> 1.7

				int j,k,l,ii;
				int energy[5];
				structure *ct;

				//convert temp to a ct file
				ii = oligo[i].length();
				ct = new structure();

				//point to location of data tables:
				ct->SetThermodynamicDataTable(data);

				ct->allocate(ii);
				//ct->numofbases=ii;

				//strcpy(ct->ctlabel[1],"");
				for (j=0;j<ii;j++) {
					if (oligo[i][j]=='A'||oligo[i][j]=='a') ct->numseq[j+1] = 1;
					else if (oligo[i][j]=='C'||oligo[i][j]=='c') ct->numseq[j+1] = 2;
					else if (oligo[i][j]=='G'||oligo[i][j]=='g') ct->numseq[j+1] = 3;
					else  ct->numseq[j+1] = 4;
				}

				//First Calculate the Bimolecular folding DG:
				energy[0] = bimolecular(ct,data);


//<<<<<<< OligoScreenCalc.cpp
		//gcvt((float (ct->energy[1]))/10,6,temp);
//		if (conversionfactor==10) sprintf(temp,"%.1f",(float (ct->GetEnergy(1)))/conversionfactor);
//		else sprintf(temp,"%.2f",(float (ct->GetEnergy(1)))/conversionfactor); //assume conversionfactor==100
//=======
//>>>>>>> 1.7

				//Now do the unimolecular folding:
				dynamic(ct,data,CNTR6,CNTR8,CNTR9);
				energy[1] = ct->GetEnergy(1);
				//efn2 (&data,ct);--no efn2 needed with current parameters
				//sortstructures (ct);
		

				//now calculate the duplex free energy:
			
//<<<<<<< OligoScreenCalc.cpp
		//if (!pObject->isRNA) {//oligo is DNA:
//		if (hybriddata!=NULL) {
//      		energy = hybriddata->init;//initiation
//			for (j=1;j<(ct->GetSequenceLength());j++) {
//          		energy = energy + hybriddata->stack[complement(j+1,ct)][ct->numseq[j+1]][complement(j,ct)][ct->numseq[j]];
//			}
//=======
				//if (!pObject->isRNA) {//oligo is DNA:
				if (hybriddata!=NULL) {
		      		energy[2] = hybriddata->init;//initiation
					for (j=1;j<(ct->GetSequenceLength());j++) {
			  		energy[2] = energy[2] + hybriddata->stack[complement(j+1,ct)][ct->numseq[j+1]][complement(j,ct)][ct->numseq[j]];
					}


				}

				else {//oligo is RNA:
		      		energy[2] = data->init;//initiation
		      		for (j=1;j<(ct->GetSequenceLength());j++) {
			  		energy[2] = energy[2] + data->stack[ct->numseq[j]]
			    		[complement(j,ct)][ct->numseq[j+1]][complement(j+1,ct)];
					}


					//add end correction for RNA-RNA duplexes
					if (ct->numseq[1]==4||ct->numseq[1]==1) {
						energy[2] = energy[2] + data->auend;
					}
					if (ct->numseq[ct->GetSequenceLength()]==4||ct->numseq[ct->GetSequenceLength()]==1) {
						energy[2] = energy[2] + data->auend;
					}

				}



				//Added 1/11/05
				//Calculate the free energy cost of opening two base pairs at 5' end of duplex and 3' end of duplex
		
				//start at 5' end of sequence submitted
		
		      		energy[3] = 0;
				for (j=1;j<(ct->GetSequenceLength())&&j<3;j++) {
			  	energy[3] = energy[3] - data->stack[ct->numseq[j]][complement(j,ct)][ct->numseq[j+1]][complement(j+1,ct)];
				}
				//check for AU/GU ends
				if (ct->numseq[1]==4||ct->numseq[1]==1) {
					energy[3] = energy[3] - data->auend;
				}
				if (ct->numseq[3]==4||ct->numseq[3]==1) {
					energy[3] = energy[3] + data->auend;
				}


				//now do 3' end of sequence submitted
				energy[4] = 0;
				for (j=ct->GetSequenceLength();j>(0)&&j>ct->GetSequenceLength()-2;j--) {
			  	energy[4] = energy[4] - data->stack[complement(j,ct)][ct->numseq[j]][complement(j-1,ct)][ct->numseq[j-1]];
				}
				//check for AU/GU ends
				if (ct->numseq[ct->GetSequenceLength()]==4||ct->numseq[ct->GetSequenceLength()]==1) {
					energy[4] = energy[4] - data->auend;
				}
				if (ct->GetSequenceLength()>2) {
					if (ct->numseq[ct->GetSequenceLength()-2]==4||ct->numseq[ct->GetSequenceLength()-2]==1) {
						energy[4] = energy[4] + data->auend;
					}
				}
		

				//write to the output file
#ifdef SMP
#pragma omp ordered
#endif
				for (int ii=0;ii<5;ii++){
					if (conversionfactor==10) out.precision(1);
					else out.precision(2); //assume conversionfactor==100
					out << fixed;
					if (ii == 0) out << oligo[i] << "\t";
			 		out << (float (energy[ii]))/conversionfactor;
					if (ii != 4) out << "\t";
					else out << "\n";
				}
//>>>>>>> 1.7

				delete ct;


//<<<<<<< OligoScreenCalc.cpp
//		else {//oligo is RNA:
//      		energy = data->init;//initiation
//      		for (j=1;j<(ct->GetSequenceLength());j++) {
//          		energy = energy + data->stack[ct->numseq[j]]
 //           		[complement(j,ct)][ct->numseq[j+1]][complement(j+1,ct)];
//=======
//>>>>>>> 1.7
			}
//<<<<<<< OligoScreenCalc.cpp


			//add end correction for RNA-RNA duplexes
//			if (ct->numseq[1]==4||ct->numseq[1]==1) {
//				energy = energy + data->auend;
//			}
//			if (ct->numseq[ct->GetSequenceLength()]==4||ct->numseq[ct->GetSequenceLength()]==1) {
//				energy = energy + data->auend;
//			}

//=======
//>>>>>>> 1.7
		}
//<<<<<<< OligoScreenCalc.cpp

		//gcvt((float (energy))/10,6,temp);
//		if (conversionfactor==10) sprintf(temp,"%.1f",(float (energy))/conversionfactor);
//		else sprintf(temp,"%.2f",(float (energy))/conversionfactor); //assume conversionfactor==100


//		out << temp<<"\t";

		//Added 1/11/05
		//Calculate the free energy cost of opening two base pairs at 5' end of duplex and 3' end of duplex
		
		//start at 5' end of sequence submitted
		
//      	energy = 0;
//		for (j=1;j<(ct->GetSequenceLength())&&j<3;j++) {
 //         	energy = energy - data->stack[ct->numseq[j]][complement(j,ct)][ct->numseq[j+1]][complement(j+1,ct)];
//		}
		//check for AU/GU ends
//		if (ct->numseq[1]==4||ct->numseq[1]==1) {
//			energy = energy - data->auend;
//		}
//		if (ct->numseq[3]==4||ct->numseq[3]==1) {
//			energy = energy + data->auend;
//		}
		
		//gcvt((float (energy))/10,6,temp);
//		if (conversionfactor==10) sprintf(temp,"%.1f",(float (energy))/conversionfactor);
///		else sprintf(temp,"%.2f",(float (energy))/conversionfactor); //assume conversionfactor==100


//		out << temp<<"\t";

		//now do 3' end of sequence submitted
		
//      	energy = 0;
//		for (j=ct->GetSequenceLength();j>(0)&&j>ct->GetSequenceLength()-2;j--) {
//          	energy = energy - data->stack[complement(j,ct)][ct->numseq[j]][complement(j-1,ct)][ct->numseq[j-1]];
//		}
		//check for AU/GU ends
//		if (ct->numseq[ct->GetSequenceLength()]==4||ct->numseq[ct->GetSequenceLength()]==1) {
//			energy = energy - data->auend;
//		}
//		if (ct->GetSequenceLength()>2) {
//			if (ct->numseq[ct->GetSequenceLength()-2]==4||ct->numseq[ct->GetSequenceLength()-2]==1) {
//				energy = energy + data->auend;
//			}
//		}
		
		//gcvt((float (energy))/10,6,temp);
//		if (conversionfactor==10) sprintf(temp,"%.1f",(float (energy))/conversionfactor);
//		else sprintf(temp,"%.2f",(float (energy))/conversionfactor); //assume conversionfactor==100


//		out << temp<<"\n";


//		delete ct;

//		in >> temp;

		
//=======
//>>>>>>> 1.7
	}

	in.close();
	out.close();
	delete[] oligo;

	return 0;//no errors
}
