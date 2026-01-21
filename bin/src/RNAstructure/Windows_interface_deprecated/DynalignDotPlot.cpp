// DynalignDotPlot.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "DynalignDotPlot.h"
#include "../src/structure.h"
#include "../src/algorithm.h"
#include "../src/dynalign.h"
#include "../src/dynalignheap.h"

using namespace std;

// CDynalignDotPlot dialog

IMPLEMENT_DYNAMIC(CDynalignDotPlot, CDialog)
CDynalignDotPlot::CDynalignDotPlot(CWnd* pParent /*=NULL*/)
	: CDialog(CDynalignDotPlot::IDD, pParent)
	, savefilename(_T(""))
	, outfilename(_T(""))
	, maxenergy (0)
	, maxpercent (0)
{
}

CDynalignDotPlot::~CDynalignDotPlot()
{
}

void CDynalignDotPlot::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_SAVEFILANAME, savefilename);
	DDX_Text(pDX, IDC_OUTFILENAME, outfilename);
	DDX_Text(pDX, IDC_MAXPERCENT, maxpercent);
	DDV_MinMaxFloat(pDX, maxpercent, 0, 100);
	DDX_Text(pDX, IDC_MAXE, maxenergy);
	DDV_MinMaxFloat(pDX, maxenergy, 0, 10000000);
}


BEGIN_MESSAGE_MAP(CDynalignDotPlot, CDialog)
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
	ON_BN_CLICKED(IDC_SAVEFILE, OnBnClickedSavefile)
	ON_BN_CLICKED(IDC_OUTFILE, OnBnClickedOutfile)
END_MESSAGE_MAP()


// CDynalignDotPlot message handlers

void CDynalignDotPlot::OnBnClickedOk()
{
	//Open the save file, make the output file
	int ir,ip,jp,kp,cp;
	short i,j,k,l;
	integersize en1;
	
	short maxsep,gap,lowest;
	//integersize ****v,****w,**w3,**w5, ****vmod,
	dynalignarray *w,*vmod;
	varray *v;
	wendarray *w3, *w5;
	integersize crit;
	structure ct1,ct2;
	datatable *data;
	dynalignheap heap;
	ofstream out;
	int modificationflag;
	data = new datatable;
	bool singleinsert,local;
	bool **allowed_alignments;
	short *lowend,*highend;

	UpdateData(TRUE);

	if(savefilename==""||outfilename=="") {
		MessageBox("Please specify names for the savefile and the output file.");
		return;
	}

	//open the save file to peek at the sizes needed to allocate arrays:
	ifstream sav(savefilename.GetBuffer(10),ios::binary);


	read(&sav, &modificationflag);

	//start with structure information
	int length1,length2;
	read(&sav,&(length1));
	read(&sav,&(length2));
	read(&sav,&maxsep);
	sav.close();
	
	if (maxsep<0) {
	//This means that the allowed alignments were constrained by the HMM alignment
	  //probs and not by the traditional M parameter.
	  //Therefore, space must be allocated for the allowed_alignment arrays.

	  allowed_alignments=new bool *[length1+1];
	  for (i=0;i<=length1;i++) allowed_alignments[i]=new bool [length2+1];

	}
	else allowed_alignments=NULL;
	//fill the low and highend arrays:
	//allocate space:
	lowend = new short [2*length1];
	highend = new short [2*length1];


	if (modificationflag==1) {
		vmod = new dynalignarray();
		

	}
	else vmod=NULL;

	v = new varray();
	w = new dynalignarray();

	

	w3 = new wendarray();
	w5 = new wendarray();

	
	opendynalignsavefile(savefilename.GetBuffer(10), &ct1, &ct2, v, w, vmod, w3, w5, data, &singleinsert,/* pair,*/ &maxsep, &gap, &lowest, &local, allowed_alignments,lowend,highend);
	
	//Now do the actual work of finding the dots

	crit = INFINITE_ENERGY;
	for (i=lowlimit(ct1.GetSequenceLength(),maxsep,ct1.GetSequenceLength(),ct2.GetSequenceLength());i</*2*maxsep+2*/highlimit(ct1.GetSequenceLength(),maxsep,ct1.GetSequenceLength(),ct2.GetSequenceLength())&&i<=ct2.GetSequenceLength();i++) {
		
			//k = i+ct1.numofbases-maxsep;
			//if (k<=ct2.numofbases) {
				if (w5->f(ct1.GetSequenceLength(),i)/*[ct1.numofbases][i]*/+gap*(ct2.GetSequenceLength()-i)/*(i-maxsep)*/<crit) crit = w5->f(ct1.GetSequenceLength(),i)/*[ct1.numofbases][i]*/+gap*(ct2.GetSequenceLength()-i);	
			//}
		
	}
	if (maxpercent>0) crit = crit - max((integersize)((float) crit*(((float) maxpercent)/100.0 )),(integersize) -maxenergy*10);
	else crit = crit - (short) maxpercent;
	for (i=1;i<=ct1.GetSequenceLength();i++) {
		for (j=i+minloop;j<=ct1.GetSequenceLength();j++) {
			for (k=max(/*i-maxsep*/lowlimit(i,maxsep,ct1.GetSequenceLength(),ct2.GetSequenceLength()),1);k<=(min(ct2.GetSequenceLength(),highlimit(i,maxsep,ct1.GetSequenceLength(),ct2.GetSequenceLength())/*i+maxsep*/));k++) {
				for (l=max(/*j-maxsep*/lowlimit(j,maxsep,ct1.GetSequenceLength(),ct2.GetSequenceLength()),k);l<=(min(ct2.GetSequenceLength(),highlimit(j,maxsep,ct1.GetSequenceLength(),ct2.GetSequenceLength())/*j+maxsep*/));l++) {
					
					//use a and b when refering to the energy arrays
					//a = k-i+maxsep;
					//b = l-j+maxsep;	

					//if ((a+ct2.numofbases-ct1.numofbases>=0)&&(a+ct2.numofbases-ct1.numofbases<=2*maxsep)) {
						if (v->f(i,j,k,l)/*v[j][i][a][b]*/+v->f(j,i+ct1.GetSequenceLength(),l,k+ct2.GetSequenceLength())/*v[i+ct1.numofbases][j-i][b][a]*/<=crit) 
							heap.push(i,j,k,l,v->f(i,j,k,l)/*v[j][i][a][b]*/+v->f(j,i+ct1.GetSequenceLength(),l,k+ct2.GetSequenceLength()));
					//}

				}
			}
		}

	}
	//make a heap out of dynalign heap:
	for (ip=1;ip<heap.size;ip++) {
		jp = ip;
		kp = jp/2;
		while ((heap.peak(jp)<heap.peak(kp))&&kp>=0){
			heap.swap(jp,kp);
			jp = jp/2;
			kp = kp/2;
		}

	}
	//sort the dynalign heap with heap sort
	for (ip = heap.size-2;ip>=0;ip--) {
		heap.swap(ip+1,0);

		jp = 0;
		cp = 1;
		while (cp<=ip) {
			if (cp!=ip) {
				if (heap.peak(cp+1)<heap.peak(cp)) cp++;

			}
			if (heap.peak(cp)<heap.peak(jp)) {
				heap.swap(cp,jp);
				jp = cp;
				cp = 2*cp;
			}
			else cp = ip + 1;
		}

	}

	//output the file:
	out.open(outfilename.GetBuffer(10));
	out << "i\tj\tk\tl\tenergy*10\n";
	for (ir=heap.size-1;ir>0;ir--) {
		heap.read(ir,&i,&j,&k,&l,&en1);
		out << i << "\t" << j << "\t" << k << "\t" << l << "\t" << en1 << "\n";
	}



	out.close();


	if (maxsep<0) {
	  
		for (i=0;i<=ct1.GetSequenceLength();i++) delete[] allowed_alignments[i];
		delete[] allowed_alignments;
	}

	if (modificationflag) delete vmod;
	delete v;
	delete w;
	delete w3;
	delete w5;

	delete lowend;
	delete highend;
		

	
	//delete[] pair[0];
	//delete[] pair[1];
	//delete[] pair;
	
	delete[] data;

	OnOK();
}

void CDynalignDotPlot::OnBnClickedSavefile()
{
	//Get the save file name:
	CFileDialog *filedialog;
	int i;
	char *name;
		
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Dynalign Save File (*.dsv)|*.dsv||");

	
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

void CDynalignDotPlot::OnBnClickedOutfile()
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
