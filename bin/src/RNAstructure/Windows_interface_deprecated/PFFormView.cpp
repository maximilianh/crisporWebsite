// PFFormView.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "PFFormView.h"
#include "doubledialog.h"
#include "pairdialog.h"
#include "singledialog.h"
#include "../src/outputconstraints.h"
#include "temp_Dialog.h"

#include "../src/algorithm.h"

#include <direct.h>
#include "../src/pfunction.h"
#include "DistanceDialog.h"
#include "readshape.h"
#include "LinearSHAPEDialog.h"


//This is the backend functions, set up to run in a seperate thread
UINT PFProc( LPVOID pParam ){    
	int i,bases;
	char *savename,buffer[7];
	double time,currenttime;


	CPFObject* pObject = (CPFObject*) pParam;
	CPFFormView* pView = (CPFFormView*) pObject->parent; 
    
	pfdatatable *data;
	data = new pfdatatable(pObject->data,scalingdefinition,pObject->Temperature);

	

	if (pObject->subfold&&!pObject->ct->intermolecular) {
		//fold each fragment of sequence that has the common 5' end, i.e. whittle way the 3' end
		bases=pObject->ct->GetSequenceLength();

		pObject->savefile[strlen(pObject->savefile)-4]='\0';
		savename=new char[strlen(pObject->savefile)+11]; //allocate enough space for sequences of length 999,999 nucs.
		time=0;
		currenttime=0;
		for (i=bases;i>minloop+1;i--) {
			time+=pow((double) i,3.0);	
		}
		

		
		structure *subct;
		for (i=bases;i>minloop+1;i--) {
			subct = new structure();
			subct->allocate(i);
			for (int ii=1;ii<=i;++ii) {
				subct->numseq[ii]=pObject->ct->numseq[ii];
				subct->nucs[ii]=pObject->ct->nucs[ii];
				subct->hnumber[ii]=pObject->ct->hnumber[ii];
			}

			//pObject->ct->GetSequenceLength()=i;
			
			
			strcpy(savename,pObject->savefile);
			strcat(savename,"_");
			sprintf(buffer,"%u",i);
			strcat(savename,buffer);
			strcat(savename,".pfs");
	
			pfunction((subct),data,0,savename);

			currenttime+=pow((double) i,3.0);
			pObject->progress->update(100.0*currenttime/time);
			delete subct;
		}
		
		//pObject->ct->numofbases=bases;

		delete[] savename;

	}
	


	else pfunction((pObject->ct),data,pObject->progress, pObject->savefile);
	
	
	
    
	::PostMessage(pView->m_hWnd,ID_FOLDDONE,0,0);
	

	return 0;   // thread completed successfully
}


// CPFFormView

IMPLEMENT_DYNCREATE(CPFFormView, CFormView)

CPFFormView::CPFFormView()
	: CFormView(CPFFormView::IDD)
	, sequencename(_T(""))
	, savename(_T(""))
{

	subfold=false;
}

CPFFormView::~CPFFormView()
{
}

void CPFFormView::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_SEQUENCENAME, sequencename);
	DDX_Text(pDX, IDC_SAVENAME, savename);
}

BEGIN_MESSAGE_MAP(CPFFormView, CFormView)
	ON_BN_CLICKED(IDC_START, OnBnClickedStart)
	ON_BN_CLICKED(IDC_SEQUENCEBUTTON, OnBnClickedSequencebutton)
	ON_BN_CLICKED(IDC_SAVEBUTTON, OnBnClickedSavebutton)
	ON_MESSAGE(ID_FOLDDONE, DoneFolding)
	ON_COMMAND(ID_FORCE_BASEPAIR, OnForceBasepair)
	ON_COMMAND(ID_FORCE_CHEMICALMODIFICATION, OnForceChemicalmodification)
	ON_COMMAND(ID_FORCE_DOUBLESTRANDED, OnForceDoublestranded)
	ON_COMMAND(ID_FORCE_FMN, OnForceFmn)
	ON_COMMAND(ID_FORCE_SINGESTRANDED, OnForceSingestranded)
	ON_COMMAND(ID_FORCE_CURRENT, OnForceCurrent)
	ON_COMMAND(ID_FORCE_RESET, OnForceReset)
	ON_COMMAND(ID_FORCE_SAVECONSTRAINTS, OnForceSaveconstraints)
	ON_COMMAND(ID_FORCE_RESTORECONSTRAINTS, OnForceRestoreconstraints)
	ON_COMMAND(ID_FORCE_FORBIDBASEPAIRS, OnForceForbidbasepairs)
	ON_COMMAND(ID_TEMPERATURE, OnTemperature)
	ON_COMMAND(ID_FORCE_MAXIMUMPAIRINGDISTANCE, OnForceMaximumpairingdistance)
	ON_COMMAND(ID_FORCE_READSHAPEREACTVITY, OnForceReadshapereactvity)
	ON_COMMAND(ID_READ_SHAPE_LINEAR, OnReadShapeLinear)
	ON_COMMAND(ID_SUBFOLD, OnSubfold)
END_MESSAGE_MAP()


// CPFFormView diagnostics

#ifdef _DEBUG
void CPFFormView::AssertValid() const
{
	CFormView::AssertValid();
}

void CPFFormView::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}
#endif //_DEBUG


// CPFFormView message handlers

void CPFFormView::OnUpdate(CView*, LPARAM, CObject*)
{

	ResizeParentToFit(FALSE);
	ResizeParentToFit();

	
	

}

void CPFFormView::OnInitialUpdate() {
	CFormView::OnInitialUpdate();

	//Below doesn't work because item does not seem to exist:	
	//CMenu* menu = GetFoldDocument()->menuframe->GetMenu( );
	//UINT ret = menu->CheckMenuItem(ID_SUBFOLD,MF_CHECKED);

}

void CPFFormView::OnBnClickedStart()
{
	
	
	//Now do the work of the partition function calculation
	UpdateData(TRUE);
	
	if(sequencename=="") {
		MessageBox("Please specify a sequence name.");
		return;
	}

	//TProgressDialog *progress;
	progress = new TProgressDialog();//parent,"Folding the RNA..."
	CRect *rect;
	rect = new CRect(10,40,210,120);
	
	if (GetFoldDocument()->ISRNA)
	progress->Create(NULL,"Folding the RNA...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);//(IDD_PROGRESSDIALOG);//
	else progress->Create(NULL,"Folding the DNA...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);
	delete rect;

	//check to see if the temperature has been changed.
	if (GetFoldDocument()->T<310||GetFoldDocument()->T>311) {

		//change the temperature from 310.15
		if (GetFoldDocument()->newtemp()==0) {
			//if newtemp returned zero, pass a warning to the user
			AfxMessageBox( "An enthalpy data file could not be found!\nTemperature of prediction will revert back to 37 degrees C.\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
				MB_OK|MB_ICONHAND);

		}

	}


	
	pfobject.parent=this;
	
	pfobject.data=&GetFoldDocument()->data;
	pfobject.savefile = new char[savename.GetLength()+1];
	strcpy(pfobject.savefile,savename.GetBuffer(10));
	pfobject.progress = progress;
	pfobject.ct=&(GetFoldDocument()->ct);
	pfobject.subfold=subfold;
	pfobject.Temperature = GetFoldDocument()->T;


	AfxBeginThread(PFProc,&pfobject);
	



}

void CPFFormView::OnBnClickedSequencebutton()
{
	//Open the standard dialog to pick a sequence
	char *tsavename;
	short int i;
	CFileDialog *filedialog;
	
	

	
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Sequence Files (*.seq)|*.seq||");

	
	filedialog->m_ofn.lpstrInitialDir=GetFoldDocument()->startpath;
	if (filedialog->DoModal()==IDOK) {
		//strcpy(m_sequencename.GetBuffer(10),(filedialog->GetPathName()).GetBuffer(0));
		sequencename=(filedialog->GetPathName()).GetBuffer(30);
		i = sequencename.GetLength();
		
		tsavename = new char[i+4];//allocate enough space so that 
														//three characters can be added 
														//to the name if necessary
		strcpy(tsavename,sequencename.GetBuffer(10));
		//count the characters to the .
		
		while(i>=0){
			
			if (tsavename[i]=='.') break;
			i--;
		}
		if (i==0) i = sequencename.GetLength();
		strcpy(tsavename+i+1,"pfs\0");
		savename=tsavename;
		
		delete[] tsavename;//fix this?
		UpdateData(FALSE);
		
		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(GetFoldDocument()->startpath,_MAX_PATH);

		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(GetFoldDocument()->startpath,_MAX_PATH);
		CString path;
		path = filedialog->GetPathName();
		i = path.GetLength();
		while(i>=0){
			
			if (path[i]=='\\') break;
			i--;
		}
		if (i>_MAX_PATH) i = _MAX_PATH;
		strncpy(GetFoldDocument()->startpath,path.GetBuffer(1),i);
		*(GetFoldDocument()->startpath + i) ='\0';

		GetFoldDocument()->ct.openseq(sequencename.GetBuffer(10));	
		

	}
	delete filedialog;
	
}

void CPFFormView::OnBnClickedSavebutton()
{
	//The user is specifying the save file name explicitly
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".pfs",savename,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"Partition Function Save Files (*.pfs)|*.pfs|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		savename=(filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);		

	}
	delete filedialog;
}

CFoldDoc *CPFFormView::GetFoldDocument() {
	
	return ((CFoldDoc*) GetDocument());	

}

LRESULT CPFFormView::DoneFolding(WPARAM wParam, LPARAM lParam) {
	delete[] pfobject.savefile;

	progress->SendMessage (WM_CLOSE);

	if (!subfold) ((CRNAstructureApp*)GetFoldDocument()->pMainFrame)->BoxPlot(savename);
	
	GetFoldDocument()->Frame->SendMessage(WM_CLOSE);

	

	
	return 0;
}

void CPFFormView::OnForceBasepair()
{
	CPairDialog *pairdialog;

	pairdialog = new CPairDialog(&(GetFoldDocument()->ct),this);

	pairdialog->Create(IDD_PAIR_DIALOG,this);
}

void CPFFormView::OnForceChemicalmodification()
{
	CSingleDialog *singledialog;

	singledialog = new CSingleDialog(&(GetFoldDocument()->ct),this,3);

	singledialog->Create(IDD_SINGLE_DIALOG,this);
}

void CPFFormView::OnForceDoublestranded()
{
	CDoubleDialog *doubledialog;

	doubledialog = new CDoubleDialog(&(GetFoldDocument()->ct),this);

	doubledialog->Create(IDD_DOUBLE_DIALOG,this);
}

void CPFFormView::OnForceFmn()
{
	CSingleDialog *singledialog;

	singledialog = new CSingleDialog(&(GetFoldDocument()->ct),this,2);

	singledialog->Create(IDD_SINGLE_DIALOG,this);
}

void CPFFormView::OnForceSingestranded()
{
	CSingleDialog *singledialog;

	singledialog = new CSingleDialog(&(GetFoldDocument()->ct),this);

	singledialog->Create(IDD_SINGLE_DIALOG,this);
}

void CPFFormView::OnForceCurrent()
{
	char message[5000],temp[5];
	int it,jt;

	strcpy (message,"");


   if (GetFoldDocument()->ct.GetNumberofDoubles()>0) {
   	strcat (message,"Bases forced double:\n");
      for (it=0;it<GetFoldDocument()->ct.GetNumberofDoubles();it++) {
      	itoa(GetFoldDocument()->ct.GetDouble(it),temp,10);
         strcat(message,"   ");
         strcat(message,temp);

         if (it%5==0) strcat(message,"\n");
		 else strcat(message,", ");
      }
	  strcat(message,"\n");
   }


   if (GetFoldDocument()->ct.GetNumberofPairs()>0) {
   	strcat(message,"Forced Base Pairs:\n");
   	for (it = 0; it<GetFoldDocument()->ct.GetNumberofPairs();it++) {

 			itoa(GetFoldDocument()->ct.GetPair5(it),temp,10);
         strcat(message,"  ");
         strcat(message,temp);
         strcat(message," - ");
         itoa(GetFoldDocument()->ct.GetPair3(it),temp,10);
         strcat(message,temp);
         if (it%5==0) strcat(message,"\n");
		 else strcat(message,", ");
      }
		strcat(message,"\n");
   }

   if (GetFoldDocument()->ct.GetNumberofSingles()>0) {
   	strcat (message,"Bases forced single:\n");
      for (it=0;it<GetFoldDocument()->ct.GetNumberofSingles();it++) {
      	itoa(GetFoldDocument()->ct.GetSingle(it),temp,10);
         strcat(message,"   ");
         strcat(message,temp);
         if (it%5==0) strcat(message,"\n");
		 else strcat(message,", ");
      }
	  strcat(message,"\n");
   }

   if (GetFoldDocument()->ct.GetNumberofGU()>0) {
   	strcat (message,"U's in GU pairs:\n");
      for (it=0;it<GetFoldDocument()->ct.GetNumberofGU();it++) {
      	itoa(GetFoldDocument()->ct.GetGUpair(it),temp,10);
         strcat(message,"   ");
         strcat(message,temp);
         if (((it+1)%5==0)) strcat(message,"\n");
		 else strcat(message,", ");
      }
	  strcat(message,"\n");
   }

   if (GetFoldDocument()->ct.GetNumberofModified()>0) {
   	strcat (message,"Bases accessible to chemical modification:\n");
      for (it=0;it<GetFoldDocument()->ct.GetNumberofModified();it++) {
      	itoa(GetFoldDocument()->ct.GetModified(it),temp,10);
         strcat(message,"   ");
         strcat(message,temp);
         if (it%5==0) strcat(message,"\n");
		 else strcat(message,", ");
      }
	  strcat(message,"\n");
   }

   if (GetFoldDocument()->ct.min_gu>0) {
	  strcat (message,"Minimum # of G-U pairs=");
      
      	itoa(GetFoldDocument()->ct.min_gu,temp,10);
         strcat(message," ");
         strcat(message,temp);
         strcat(message,"\n");
		 
   }

   if (GetFoldDocument()->ct.neighbors[0][0]>0) {
	   strcat (message,"Paired neighbors:");
      
      for (it=0;it<100&&GetFoldDocument()->ct.neighbors[it][0]>0;it++) {
		  for (jt=0;jt<25&&GetFoldDocument()->ct.neighbors[it][jt]>0;jt++) {
			strcat(message,tobase(GetFoldDocument()->ct.neighbors[it][jt]));
		  }
		  strcat(message,"\n");
      }
	  strcat(message,"\n");
		 
   }

   if (GetFoldDocument()->ct.min_gu>0) {
	  strcat (message,"Minimum # of G-U pairs=");
      
      	itoa(GetFoldDocument()->ct.min_gu,temp,10);
         strcat(message," ");
         strcat(message,temp);
         strcat(message,"\n");
		 
   }

   if (GetFoldDocument()->ct.min_g_or_u>0) {
	  strcat (message,"Minimum # of G's and U's in pairs=");
      
      	itoa(GetFoldDocument()->ct.min_g_or_u,temp,10);
         strcat(message," ");
         strcat(message,temp);
         strcat(message,"\n");
		 
   }

   if (GetFoldDocument()->ct.GetNumberofForbiddenPairs()>0) {
		strcat (message,"Prohibited Base Pairs=");
		for (it = 0; it<GetFoldDocument()->ct.GetNumberofForbiddenPairs();it++) {

 			itoa(GetFoldDocument()->ct.GetForbiddenPair5(it),temp,10);
			strcat(message,"  ");
			strcat(message,temp);
			strcat(message," - ");
			itoa(GetFoldDocument()->ct.GetForbiddenPair3(it),temp,10);
			strcat(message,temp);
			if (it%5==0) strcat(message,"\n");
			else strcat(message,", ");
		}
		strcat(message,"\n");
   }

   

	if (GetFoldDocument()->ct.GetNumberofDoubles()==0&&GetFoldDocument()->ct.GetNumberofPairs()==0&&
		GetFoldDocument()->ct.GetNumberofSingles()==0&&GetFoldDocument()->ct.GetNumberofGU()==0&&
		GetFoldDocument()->ct.GetNumberofModified()==0&&GetFoldDocument()->ct.min_gu==0&&
		GetFoldDocument()->ct.min_g_or_u==0&&GetFoldDocument()->ct.neighbors[0][0]==0
		&&GetFoldDocument()->ct.GetNumberofForbiddenPairs()==0) {

		strcpy(message,"There are no folding constraints.");
	}

   AfxMessageBox( message, 
			MB_OK|MB_ICONINFORMATION   );
}

void CPFFormView::OnForceReset()
{
	if (AfxMessageBox( "This will erase all constraints.\nContinue?", 
		MB_OKCANCEL   |MB_ICONQUESTION   )==IDOK) {
		//reset all constraints:
		//GetFoldDocument()->ct.ndbl = 0;
		//GetFoldDocument()->ct.npair = 0;
		//GetFoldDocument()->ct.nnopair = 0;
		//GetFoldDocument()->ct.nmod = 0;
		//GetFoldDocument()->ct.ngu=0;
		//GetFoldDocument()->ct.nforbid=0;
			GetFoldDocument()->ct.RemoveConstraints();
	}
}

void CPFFormView::OnForceSaveconstraints()
{
	//The user is specifying the CT file name explicitly
	CFileDialog *filedialog;


	filedialog = new CFileDialog(FALSE,".con",NULL,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"Constraint Files (*.con)|*.con|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {


		outputconstraints(filedialog->GetPathName().GetBuffer(0),&(GetFoldDocument()->ct));

	

	}
	delete filedialog;
}

void CPFFormView::OnForceRestoreconstraints()
{
	CFileDialog *filedialog;
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Constraint Files (*.con)|*.con|All Files|*.*||");

	
	
	if (filedialog->DoModal()==IDOK) {
		readconstraints(filedialog->GetPathName().GetBuffer(0),&(GetFoldDocument()->ct));
		

	}
	delete filedialog;
}

void CPFFormView::OnForceForbidbasepairs()
{
	//User wants to forbid a basepair
	CPairDialog *pairdialog;

	pairdialog = new CPairDialog(&(GetFoldDocument()->ct),this,true);

	pairdialog->Create(IDD_PAIR_DIALOG,this);

}

void CPFFormView::OnTemperature()
{
	//Allow the user to specify a new temperature.
	CTemp_Dialog *temp;
	temp=new CTemp_Dialog(&(GetFoldDocument()->T));

	temp->DoModal();
	delete temp;
	


}

void CPFFormView::OnForceMaximumpairingdistance()
{
	//Use has chosen to set/change the maximum pairing distance:
	CDistanceDialog *dist;
	dist = new CDistanceDialog(&GetFoldDocument()->ct,this);
	dist->DoModal();
	delete dist;
}

void CPFFormView::OnForceReadshapereactvity()
{
	//Use has chosen to read SHAPE reactivity data:
	CReadSHAPE *readshape;
	readshape = new CReadSHAPE(&GetFoldDocument()->ct,this);
	readshape->DoModal();
	delete readshape;
}

void CPFFormView::OnReadShapeLinear()
{
	//Use has chosen to read SHAPE reactivity data:
	CLinearSHAPEDialog *readshape;
	readshape = new CLinearSHAPEDialog(&GetFoldDocument()->ct,this);
	readshape->DoModal();
	delete readshape;
}

void CPFFormView::OnSubfold()
{
	//User has chosen the SubFold menu option	
	CMenu* menu = GetFoldDocument()->menuframe->GetMenu( );	

	subfold = !subfold;

	if (subfold)  menu->CheckMenuItem(ID_SUBFOLD,MF_CHECKED);
	else  menu->CheckMenuItem(ID_SUBFOLD,MF_UNCHECKED); 
}
