// DotPlot.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "DotPlot.h"
#include "../src/algorithm.h"

using namespace std;

// CDotPlot dialog

IMPLEMENT_DYNAMIC(CDotPlot, CDialog)
CDotPlot::CDotPlot(CWnd* pParent /*=NULL*/)
	: CDialog(CDotPlot::IDD, pParent)
	, savefilename(_T(""))
	, outfilename(_T(""))
	, maxenergy (0)
	, maxpercent (0)
{
}

CDotPlot::~CDotPlot()
{
}

void CDotPlot::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_SAVEFILANAME, savefilename);
	DDX_Text(pDX, IDC_OUTFILENAME, outfilename);
	DDX_Text(pDX, IDC_MAXPERCENT, maxpercent);
	DDV_MinMaxFloat(pDX, maxpercent, 0, 100);
	DDX_Text(pDX, IDC_MAXE, maxenergy);
	DDV_MinMaxFloat(pDX, maxenergy, 0, 10000000);
}


BEGIN_MESSAGE_MAP(CDotPlot, CDialog)
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
	ON_BN_CLICKED(IDC_SAVEFILE, OnBnClickedSavefile)
	ON_BN_CLICKED(IDC_OUTFILE, OnBnClickedOutfile)
END_MESSAGE_MAP()


// CDotPlot message handlers

void CDotPlot::OnBnClickedOk()
{
	int i,j,sort;

	register int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	{0,1,0,1,0,0},{0,0,0,0,0,0}};
		
	DynProgArray<integersize> *w2,*wmb2;
	integersize *w5,*w3;
	int vmin;
	bool *lfce,*mod;
	datatable *data;
	int *heapi,*heapj;
	integersize *energy;
	int num,crit;
	int cur,c,cntr;
	structure ct;
	int q,up,ir;



	//OK calculate the dots:
	UpdateData(TRUE);

	if(savefilename==""||outfilename=="") {
		MessageBox("Please specify names for the savefile and the output file.");
		return;
	}

	

	data = new datatable;

	//peek at the length of the sequence and whether the folding is intermolecular to allocate arrays:
	ifstream sav(savefilename.GetBuffer(10),ios::binary);
	
	//Read the save file version information.
	short vers;
	read(&sav,&vers);

	//check the version
	if (vers!=safiversion) {
		//Wrong version!
		sav.close();

		MessageBox("Incorrect save file version. Please refold the sequence first with this version of RNAstructure.");

		return;
	}

	int sequencelength;
	read(&sav,&(sequencelength));
	read(&sav,&(ct.intermolecular));
	sav.close();
	//read(&sav,&(ct.numofbases));
	//read(&sav,&(ct.intermolecular));
	//sav.close();


	//allocate everything
	ct.allocate(sequencelength);
	
	DynProgArray<integersize> w(ct.GetSequenceLength());
	DynProgArray<integersize> v(ct.GetSequenceLength());
	DynProgArray<integersize> wmb(ct.GetSequenceLength());	
	forceclass fce(ct.GetSequenceLength());


	lfce = new bool [2*ct.GetSequenceLength()+1];
	mod = new bool [2*ct.GetSequenceLength()+1];
	
	w5 = new integersize [ct.GetSequenceLength()+1];
	w3 = new integersize [ct.GetSequenceLength()+2];

	if (ct.intermolecular) {
		w2 = new DynProgArray(ct.GetSequenceLength());
		wmb2 = new DynProgArray(ct.GetSequenceLength());
		
		for (i=0;i<3;i++) read(&sav,&(ct.inter[i]));
		
	}
	else {
		w2 = NULL;
		wmb2 = NULL;
	}
	
	readsav(savefilename.GetBuffer(10), &ct, w2, wmb2, w5, w3, lfce, mod, data,
			 &v, &w, &wmb, &fce, &vmin);

	//find the dots here:
	sort = 100000;
	energy = new integersize [sort+1];
	heapi = new int [sort+1];
	heapj = new int [sort+1];

	crit= min((short) (abs(vmin)*(float (maxpercent)/100.0)),(integersize) maxenergy*10);

	crit = crit + vmin;

	num = 0;
	i = 1;
	j = 2;
	while (i<(ct.GetSequenceLength())) {
		if (num==sort) {
			//allocate more space for the heap
   			delete[] heapi;
			delete[] heapj;
			delete[] energy;
			sort = 10*sort;
			heapi = new int [sort+1];
   			heapj = new int [sort+1];
   			energy = new integersize [sort+1];
		i = 1;
		j = 2;
		num = 0;

		}



		if ((v.f(i,j)+v.f(j,i+ct.GetSequenceLength()))<=crit) {

   			num++;
   			heapi[num]=i;
   			heapj[num]=j;
			energy[num] = (v.f(i,j)+v.f(j,i+ct.GetSequenceLength()));
   			j++;
   			if (j>ct.GetSequenceLength()) {
   				i++;
      			j=i+1;
   			}
		}
		else if (mod[i]||mod[j]) {
			//add i-j pair to heap if it is adjacent to unpaired nucs in one direction
			if (v.f(i,j)<INFINITE_ENERGY) {
				cur = v.f(i,j)+v.f(j+1,i+ct.GetSequenceLength()-1)+erg1(j,i+ct.GetSequenceLength(),j+1,i+ct.GetSequenceLength()-1,&ct,data);
				if (cur<=crit) {
					num++;
					heapi[num]=i;
					heapj[num]=j;
					energy[num] = crit;
   					j++;
   					if (j>ct.GetSequenceLength()) {
   						i++;
      					j=i+1;
   					}

				}
				else {
					j++;
					if (j>ct.GetSequenceLength()) {
      					i++;
						j=i+1;
					}

				}

			}
			else if (v.f(j,i+ct.GetSequenceLength())==INFINITE_ENERGY) {
				cur = v.f(i+1,j-1)+v.f(j,i+ct.GetSequenceLength())+erg1(i,j,i+1,j-1,&ct,data);
				if (cur<=crit) {
					num++;
					heapi[num]=i;
					heapj[num]=j;
					energy[num] = crit;
   					j++;
   					if (j>ct.GetSequenceLength()) {
   						i++;
      					j=i+1;
   					}

				}
				else {
					j++;
					if (j>ct.GetSequenceLength()) {
      					i++;
						j=i+1;
					}

				}

			}
			else {
				j++;
				if (j>ct.GetSequenceLength()) {
      				i++;
					j=i+1;
				}

			}

		}
		else {
   			j++;
			if (j>ct.GetSequenceLength()) {
      			i++;
				j=i+1;
			}
		}
	}



	//sort the base pair list:


	///////////////////////////////////

	//make a heap:

	
	for (q=2;q<=num;q++) {
		cur = q;
		up = cur/2;
		while ((energy[cur]<energy[up])&&up>=1) {
			swap(&heapi[cur],&heapi[up]);
			swap(&heapj[cur],&heapj[up]);
			swap(&energy[cur],&energy[up]);
			cur = cur/2;
			up = up/2;
		}
	}

	//sort the heap:

	for (ir=num-1;ir>=1;ir--) {
		swap(&heapi[ir+1],&heapi[1]);
		swap(&heapj[ir+1],&heapj[1]);
		swap(&energy[ir+1],&energy[1]);

		up =1 ;
		c = 2;
		while (c<=ir) {
			if (c!=ir) {
				if (energy[c+1]<energy[c]) c++;
			}
			if (energy[c]<energy[up]) {
				swap(&heapi[c],&heapi[up]);
				swap(&heapj[c],&heapj[up]);
				swap(&energy[c],&energy[up]);
				up=c;
				c=2*c;
			}
			else c = ir+1;
		}
	}

	ofstream out(outfilename.GetBuffer(10));
	out << "i\tj\tenergy\n";
	for(cntr = num;cntr>0;cntr--) {
		
		out << heapi[cntr]<<"\t"<<heapj[cntr]<<"\t"<<energy[cntr]<<"\n";

	}
	out.close();


	
	delete[] lfce;
	delete[] mod;
	



	delete[] w5;
	delete[] w3;


	if (ct.intermolecular) {
		delete w2;
		delete wmb2;
	}

	delete data;




	OnOK();
}

void CDotPlot::OnBnClickedSavefile()
{
		//Get the save file name:
	CFileDialog *filedialog;
	int i;
	char *name;
		
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Dynalign Save File (*.sav)|*.sav||");

	
	//filedialog->m_ofn.lpstrInitialDir=GetFoldDocument()->startpath;
	if (filedialog->DoModal()==IDOK) {
		
		savefilename=(filedialog->GetPathName()).GetBuffer(30);
		i = savefilename.GetLength();
		
		name = new char[i+4];//allocate enough space so that 
														//three characters can be added 
														//to the name if necessary
		strcpy(name,savefilename.GetBuffer(10));
		//count the characters to the .
		
		while(i>=0){
			
			if (name[i]=='.') break;
			i--;
		}
		if (i==0) i = savefilename.GetLength();
		strcpy(name+i+1,"out\0");
		outfilename=name;
		
		delete[] name;
		UpdateData(FALSE);
		
		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(GetFoldDocument()->startpath,_MAX_PATH);



		
	}
	delete filedialog;
}

void CDotPlot::OnBnClickedOutfile()
{
		// The user want to override the output file name:
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".out",outfilename,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"Out Files (*.out)|*.out|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		outfilename=(filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);		

	}
	delete filedialog;
}
