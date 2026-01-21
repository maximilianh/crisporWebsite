// FoldView.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "FoldView.h"
#include "globals.h"
#include "advanced.h"
#include "doubledialog.h"
#include "pairdialog.h"
#include "singledialog.h"
#include "../src/outputconstraints.h"
#include "nmrdialog.h"
#include "temp_Dialog.h"
#include "internaldialog.h"
#include "readshape.h"
#include "DistanceDialog.h"
#include "LinearSHAPEDialog.h"


#include <direct.h>



#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif






//This is the backend functions, set up to run in a seperate thread
UINT FoldProc( LPVOID pParam ){    
	CFoldObject* pObject = (CFoldObject*) pParam;
	CFoldView* pView = (CFoldView*) pObject->parent; 
	int i,bases;
	char *ctname,buffer[7];
	double time,currenttime;
	
    
	

	if (pObject->subfold) {
		//fold each fragment of sequence that has the common 5' end, i.e. whittle way the 3' end
		bases=pObject->ct->GetSequenceLength();

		pObject->ctoutfile[strlen(pObject->ctoutfile)-3]='\0';
		ctname=new char[strlen(pObject->ctoutfile)+10]; //allocate enough space for sequences of length 999,999 nucs.
		time=0;
		currenttime=0;
		for (i=bases;i>minloop+1;i--) {
			time+=pow((double) i,3.0);	
		}
		

		
		structure *shortct;
		for (i=bases;i>minloop+1;i--) {
			

			shortct=new structure();
			shortct->allocate(i);
			for (int ii=1;ii<=i;++ii) {
				shortct->numseq[ii]=pObject->ct->numseq[ii];
				shortct->nucs[ii]=pObject->ct->nucs[ii];
				shortct->hnumber[ii]=pObject->ct->hnumber[ii];
			}

			//pObject->ct->numofbases=i;
			dynamic(shortct,pObject->data,pObject->dynnumber,pObject->dynpercent,
				pObject->dynwindow,0,false,0,pObject->maximuminternalloopsize);
			
			strcpy(ctname,pObject->ctoutfile);
			strcat(ctname,"_");
			sprintf(buffer,"%u",i);
			strcat(ctname,buffer);
			strcat(ctname,".ct");
	
			shortct->ctout(ctname);

			currenttime+=pow((double) i,3.0);
			pObject->progress->update(100.0*currenttime/time);
			delete shortct;
		}
		
		//pObject->ct->numofbases=bases;
		delete[] ctname;
	}
	else {
		dynamic(pObject->ct,pObject->data,pObject->dynnumber,pObject->dynpercent,
			pObject->dynwindow,pObject->progress,false,pObject->savefile,pObject->maximuminternalloopsize);

	
		pObject->ct->ctout(pObject->ctoutfile);

	}
      
	::PostMessage(pView->m_hWnd,ID_FOLDDONE,0,0);
	

	return 0;   // thread completed successfully
}

//////////////////////////////////////////////////////////////////////////
void currentconstraints(structure *ct) {
	char message[5000],temp[5];
	int it,jt,i;

	strcpy (message,"");


   if (ct->GetNumberofDoubles()>0) {
   	strcat (message,"Bases forced double:\n");
      for (it=0;it<ct->GetNumberofDoubles();it++) {
      	itoa(ct->GetDouble(it),temp,10);
         strcat(message,"   ");
         strcat(message,temp);

         if (it%5==0) strcat(message,"\n");
		 else strcat(message,", ");
      }
	  strcat(message,"\n");
   }


   if (ct->GetNumberofPairs()>0) {
   	strcat(message,"Forced Base Pairs:\n");
   	for (it = 0; it<ct->GetNumberofPairs();it++) {

 		itoa(ct->GetPair5(it),temp,10);
         strcat(message,"  ");
         strcat(message,temp);
         strcat(message," - ");
         itoa(ct->GetPair3(it),temp,10);
         strcat(message,temp);
         if (it%5==0) strcat(message,"\n");
		 else strcat(message,", ");
      }
		strcat(message,"\n");
   }

   if (ct->GetNumberofSingles()>0) {
   	strcat (message,"Bases forced single:\n");
      for (it=0;it<ct->GetNumberofSingles();it++) {
      	itoa(ct->GetSingle(it),temp,10);
         strcat(message,"   ");
         strcat(message,temp);
         if (it%5==0) strcat(message,"\n");
		 else strcat(message,", ");
      }
	  strcat(message,"\n");
   }

   if (ct->GetNumberofGU()>0) {
   	strcat (message,"U's in GU pairs:\n");
      for (it=0;it<ct->GetNumberofGU();it++) {
      	itoa(ct->GetGUpair(it),temp,10);
         strcat(message,"   ");
         strcat(message,temp);
         if (((it+1)%5==0)) strcat(message,"\n");
		 else strcat(message,", ");
      }
	  strcat(message,"\n");
   }

   if (ct->GetNumberofModified()>0) {
   	strcat (message,"Bases accessible to chemical modification:\n");
      for (it=0;it<ct->GetNumberofModified();it++) {
      	itoa(ct->GetModified(it),temp,10);
         strcat(message,"   ");
         strcat(message,temp);
         if (it%5==0) strcat(message,"\n");
		 else strcat(message,", ");
      }
	  strcat(message,"\n");
   }

   

   if (ct->neighbors[0][0]>0) {
	   strcat (message,"Paired neighbors:");
      
      for (it=0;it<100&&ct->neighbors[it][0]>0;it++) {
		  for (jt=0;jt<25&&ct->neighbors[it][jt]>0;jt++) {
			strcat(message,tobase(ct->neighbors[it][jt]));
		  }
		  strcat(message,"\n");
      }
	  strcat(message,"\n");
		 
   }

   if (ct->min_gu>0) {
	  strcat (message,"Minimum # of G-U pairs=");
      
      	itoa(ct->min_gu,temp,10);
         strcat(message," ");
         strcat(message,temp);
         strcat(message,"\n");
		 
   }

   if (ct->min_g_or_u>0) {
	  strcat (message,"Minimum # of G's and U's in pairs=");
      
      	itoa(ct->min_g_or_u,temp,10);
         strcat(message," ");
         strcat(message,temp);
         strcat(message,"\n");
		 
   }

   if (ct->GetNumberofForbiddenPairs()>0) {
		strcat (message,"Prohibited Base Pairs=");
		for (it = 0; it<ct->GetNumberofForbiddenPairs();it++) {

 			itoa(ct->GetForbiddenPair5(it),temp,10);
			strcat(message,"  ");
			strcat(message,temp);
			strcat(message," - ");
			itoa(ct->GetForbiddenPair3(it),temp,10);
			strcat(message,temp);
			if (it%5==0) strcat(message,"\n");
			else strcat(message,", ");
		}
		strcat(message,"\n");
   }

   for (i=0;i<ct->nregion;i++) {

	    strcat(message,"Region from ");
		itoa(ct->start[i],temp,10);
		strcat(message,temp);
		strcat(message," to ");
		itoa(ct->stop[i],temp,10);
		strcat(message,temp);
		strcat(message,":\n");
		if (ct->rneighbors[i][0][0]>0) {
			strcat (message,"\tPaired neighbors:");
      
			for (it=0;it<100&&ct->rneighbors[i][it][0]>0;it++) {
				for (jt=0;jt<25&&ct->rneighbors[i][it][jt]>0;jt++) {
					strcat(message,tobase(ct->rneighbors[i][it][jt]));
				}
				strcat(message,"\n");
			 }
			 strcat(message,"\n");
		 
			}

		if (ct->rmin_gu[i]>0) {
			strcat (message,"\tMinimum # of G-U pairs=");
      
      		itoa(ct->rmin_gu[i],temp,10);
			strcat(message," ");
			strcat(message,temp);
			strcat(message,"\n");
		 
		}

		if (ct->rmin_g_or_u[i]>0) {
			strcat (message,"\tMinimum # of G's and U's in pairs=");
      
      		itoa(ct->rmin_g_or_u[i],temp,10);
			strcat(message," ");
			strcat(message,temp);
			strcat(message,"\n");
		 
		}

   }

   for (i=0;i<ct->nmicroarray;i++) {
		strcat(message,"Region from ");
		itoa(ct->microstart[i],temp,10);
		strcat(message,temp);
		strcat(message," to ");
		itoa(ct->microstop[i],temp,10);
		strcat(message,temp);
		strcat(message," has a minimum of ");
		itoa(ct->microunpair[i],temp,10);
		strcat(message,temp);
		strcat(message," unpaired nucs.\n");


   }

   

	if (ct->GetNumberofDoubles()==0&&ct->GetNumberofPairs()==0&&
		ct->GetNumberofSingles()==0&&ct->GetNumberofGU()==0&&
		ct->GetNumberofModified()==0&&ct->min_gu==0&&
		ct->min_g_or_u==0&&ct->nneighbors==0
		&&ct->GetNumberofForbiddenPairs()==0&&ct->nregion==0&&ct->nmicroarray==0) {

		strcpy(message,"There are no folding constraints.");
	}

   AfxMessageBox( message, 
			MB_OK|MB_ICONINFORMATION   );
}

void resetconstraints(structure *ct) {
		ct->RemoveConstraints();
		//ct->ndbl = 0;
		//ct->npair = 0;
		//ct->nnopair = 0;
		//ct->nmod = 0;
		//ct->ngu=0;
		//ct->nforbid=0;
		ct->min_gu=0;
		ct->min_g_or_u=0;
		ct->nneighbors=0;
		ct->nregion=0;
		ct->nmicroarray=0;
		
}


/////////////////////////////////////////////////////////////////////////////
// CFoldView

IMPLEMENT_DYNCREATE(CFoldView, CFormView)

CFoldView::CFoldView()
	: CFormView(CFoldView::IDD)
{

	
	//{{AFX_DATA_INIT(CFoldView)
	m_ctname = _T("");
	m_number = 20;
	m_percent = 10;
	m_save = FALSE;
	m_window = 0;
	m_sequencename = _T("");
	//}}AFX_DATA_INIT

	started=false;
	subfold=false;

	

}

void CFoldView::OnInitialUpdate() {
	char *ctname;
	int i;

	//turn off subfold (below does not work because menu does not yet exist...)	
	//CMenu* menu = GetFoldDocument()->menuframe->GetMenu( );
	//menu->CheckMenuItem(ID_SUBFOLD,MF_UNCHECKED);

	if (GetFoldDocument()->filenamedefined) {

		m_sequencename = GetFoldDocument()->filename;

		i = m_sequencename.GetLength();
		
		ctname = new char[i+4];//allocate enough space so that 
														//three characters can be added 
														//to the name if necessary
		strcpy(ctname,m_sequencename.GetBuffer(10));
		//count the characters to the .
		
		while(i>=0){
			
			if (ctname[i]=='.') break;
			i--;
		}
		if (i==0) i = m_sequencename.GetLength();
		strcpy(ctname+i+1,"ct\0");
		m_ctname=ctname;
		
		delete[] ctname;//fix this?

		NewSequence();


	}

	if (GetFoldDocument()->savefile!=NULL) {
		if (*GetFoldDocument()->savefile) {
			m_save = TRUE;
			UpdateData(FALSE);

		}
	}
	ResizeParentToFit();
	

}

CFoldView::~CFoldView()
{
}

void CFoldView::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CFoldView)
	DDX_Text(pDX, IDC_CTNAME, m_ctname);
	DDX_Text(pDX, IDC_NUMBER, m_number);
	DDX_Text(pDX, IDC_PERCENT, m_percent);
	DDV_MinMaxInt(pDX, m_percent, 0, 99);
	DDX_Check(pDX, IDC_SAVECHECK, m_save);
	DDX_Text(pDX, IDC_WINDOW, m_window);
	DDV_MinMaxInt(pDX, m_window, 0, 16960);
	DDX_Text(pDX, IDC_SEQUENCENAME, m_sequencename);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CFoldView, CFormView)
	ON_MESSAGE(ID_FOLDDONE, DoneFolding)
	//{{AFX_MSG_MAP(CFoldView)
	ON_BN_CLICKED(IDC_SEQUENCEBUTTON, OnSequencebutton)
	ON_BN_CLICKED(IDC_CTBUTTON, OnCtbutton)
	ON_BN_CLICKED(IDC_START, OnStart)
	ON_COMMAND(ID_FORCE_BASEPAIR, OnForceBasepair)
	ON_COMMAND(ID_FORCE_CURRENT, OnForceCurrent)
	ON_COMMAND(ID_FORCE_DOUBLESTRANDED, OnForceDoublestranded)
	ON_COMMAND(ID_FORCE_RESET, OnForceReset)
	ON_COMMAND(ID_FORCE_RESTORECONSTRAINTS, OnForceRestoreconstraints)
	ON_COMMAND(ID_FORCE_SAVECONSTRAINTS, OnForceSaveconstraints)
	ON_COMMAND(ID_FORCE_SINGESTRANDED, OnForceSingestranded)
	//ON_COMMAND(ID_ADVANCED_SUBOPTIMALSTRUCTURES, OnAdvancedSuboptimalstructures)
	ON_COMMAND(ID_FORCE_FMN, OnForceFmn)
	ON_COMMAND(ID_FORCE_CHEMICALMODIFICATION,OnForceChemical)
	ON_EN_CHANGE(IDC_PERCENT, OnChangeParameter)
	ON_EN_CHANGE(IDC_NUMBER, OnChangeParameter)
	ON_EN_CHANGE(IDC_WINDOW, OnChangeParameter)
	//}}AFX_MSG_MAP
	ON_COMMAND(ID_FORCE_NMRCONSTRAINTS, OnForceNmrconstraints)
	ON_COMMAND(ID_FORCE_FORBIDBASEPAIRS, OnForceForbidbasepairs169)
	ON_COMMAND(ID_TEMPERATURE, OnTemperature)
	ON_COMMAND(ID_MAXIMUMLOOP, OnMaximumloop)
	ON_COMMAND(ID_FORCE_READSHAPEREACTVITY, OnForceReadshapereactvity)
	ON_COMMAND(ID_FORCE_MAXIMUMPAIRINGDISTANCE, OnForceMaximumpairingdistance)
	ON_COMMAND(ID_READ_SHAPE_LINEAR, OnReadShapeLinear)
	ON_COMMAND(ID_SUBFOLD, OnSubfold)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CFoldView diagnostics

#ifdef _DEBUG
void CFoldView::AssertValid() const
{
	CFormView::AssertValid();
}

void CFoldView::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CFoldView message handlers


void CFoldView::OnUpdate(CView*, LPARAM, CObject*)
{

	ResizeParentToFit(FALSE);
	ResizeParentToFit();
	

}


void CFoldView::OnSequencebutton() 
{
	//Open the standard dialog to pick a sequence
	char *ctname;
	short int i;
	CFileDialog *filedialog;
	
	

	
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Sequence Files (*.seq)|*.seq||");

	
	filedialog->m_ofn.lpstrInitialDir=GetFoldDocument()->startpath;
	if (filedialog->DoModal()==IDOK) {
		//strcpy(m_sequencename.GetBuffer(10),(filedialog->GetPathName()).GetBuffer(0));
		m_sequencename=(filedialog->GetPathName()).GetBuffer(30);
		i = m_sequencename.GetLength();
		
		ctname = new char[i+4];//allocate enough space so that 
														//three characters can be added 
														//to the name if necessary
		strcpy(ctname,m_sequencename.GetBuffer(10));
		//count the characters to the .
		
		while(i>=0){
			
			if (ctname[i]=='.') break;
			i--;
		}
		if (i==0) i = m_sequencename.GetLength();
		strcpy(ctname+i+1,"ct\0");
		m_ctname=ctname;
		
		delete[] ctname;//fix this?
		UpdateData(FALSE);
		
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


		NewSequence();
		

	}
	delete filedialog;

	
}


//call this function when a new sequence is chosen to set the default window size;
void CFoldView::NewSequence() 
{

	
	
	//load the ct file:
	if((GetFoldDocument()->ct.openseq(m_sequencename.GetBuffer(10)))==0) {
   	//errmsg(5,5);
   }
   if ((GetFoldDocument()->ct.GetSequenceLength())>1200) {
   			m_window=20;
   }
   else if ((GetFoldDocument()->ct.GetSequenceLength())>800) {
   			m_window=15;
            
   }
   else if ((GetFoldDocument()->ct.GetSequenceLength())>500) {
   			m_window=11;
            
   }
   else if ((GetFoldDocument()->ct.GetSequenceLength())>300) {
   			m_window=7;
            
   }
   else if ((GetFoldDocument()->ct.GetSequenceLength())>120) {
   			m_window=5;
            
   }
   else if ((GetFoldDocument()->ct.GetSequenceLength())>50) {
   			m_window=3;
            
   }
	else m_window=2;



	UpdateData(FALSE);
}


void CFoldView::OnCtbutton() 
{
	//The user is specifying the CT file name explicitly
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".ct",m_ctname,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"CT Files (*.ct)|*.ct|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		m_ctname=(filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);		

	}
	delete filedialog;
	
}

void CFoldView::OnStart() 
{
	int i;

	//a trap to make sure the calculation is started more than once
	if (started) return;
	started=true;

	if(m_sequencename=="") {
		MessageBox("Please specify a sequence name.");
		return;
	}
	UpdateData(TRUE);



    //if (dynnumber<m_number) dynnumber = m_number;
    //if (dynpercent<m_percent) dynpercent = m_percent;
    //if (dynwindow>m_window) dynwindow = m_window;
	//check to see if the temperature has been changed.
	if (GetFoldDocument()->T<310||GetFoldDocument()->T>311) {

		//change the temperature from 310.15
		if (GetFoldDocument()->newtemp()==0) {
			//if newtemp returned zero, pass a warning to the user
			AfxMessageBox( "An enthalpy data file could not be found!\nTemperature of prediction will revert back to 37 degrees C.\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
				MB_OK|MB_ICONHAND);

		}

	}

	//TProgressDialog *progress;
	progress = new TProgressDialog();//parent,"Folding the RNA..."
	CRect *rect;
	rect = new CRect(10,40,210,120);
	
	if (GetFoldDocument()->ISRNA)
	progress->Create(NULL,"Folding the RNA...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);//(IDD_PROGRESSDIALOG);//
	else progress->Create(NULL,"Folding the DNA...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);
	delete rect;


	
	
	foldobject.parent=this;
	foldobject.ct=&(GetFoldDocument()->ct);
	foldobject.data=&GetFoldDocument()->data;
	foldobject.dynpercent = m_percent;
	foldobject.dynwindow = m_window;
	foldobject.progress = progress;
	foldobject.ctoutfile = m_ctname.GetBuffer(10);
	foldobject.dynnumber = m_number ;
	foldobject.maximuminternalloopsize=GetFoldDocument()->maximuminternalloopsize;

	if (m_save) {
		i = m_ctname.GetLength();
		
		foldobject.savefile = new char[i+4];
		strcpy(foldobject.savefile,m_ctname.GetBuffer(10));
		//count the characters to the .
		
		while(i>=0){
			
			if (foldobject.savefile[i]=='.') break;
			i--;
		}
		if (i==0) i = strlen(foldobject.savefile);
		strcpy(foldobject.savefile+i+1,"sav\0");
		*GetFoldDocument()->savefile = true; 

	}
	else {
		foldobject.savefile = 0;
		*GetFoldDocument()->savefile = false; 
	}

	
	foldobject.subfold=subfold;

	
	
	AfxBeginThread(FoldProc,&foldobject);
	

	


	

	
}

LRESULT CFoldView::DoneFolding(WPARAM wParam, LPARAM lParam) {
	CDialog* finished;
	//CMainFrame* frame = (CMainFrame*) Parent;

	if (m_save) delete foldobject.savefile;

	if (!subfold) {
		//offer to display the predicted structures if this is not subfolding
		finished = new CDialog(IDD_FINISHED,this);

		if(finished->DoModal()==IDOK) ((CRNAstructureApp*) GetFoldDocument()->pMainFrame)->Draw(m_ctname.GetBuffer(10));
	

		delete finished;
	}

	

	progress->SendMessage (WM_CLOSE);
	
	GetFoldDocument()->Frame->SendMessage(WM_CLOSE);
	
	return 0;


}


CFoldDoc *CFoldView::GetFoldDocument() {
	
	return ((CFoldDoc*) GetDocument());	

}

void CFoldView::OnForceBasepair() 
{
	CPairDialog *pairdialog;

	pairdialog = new CPairDialog(&(GetFoldDocument()->ct),this);

	pairdialog->Create(IDD_PAIR_DIALOG,this);
	
}

void CFoldView::OnForceCurrent() 
{
	
	
	currentconstraints(&(GetFoldDocument()->ct));
	
	
}

void CFoldView::OnForceDoublestranded() 
{
	CDoubleDialog *doubledialog;

	doubledialog = new CDoubleDialog(&(GetFoldDocument()->ct),this);

	doubledialog->Create(IDD_DOUBLE_DIALOG,this);
	
}

void CFoldView::OnForceReset() 
{
	if (AfxMessageBox( "This will erase all constraints.\nContinue?", 
		MB_OKCANCEL   |MB_ICONQUESTION   )==IDOK) {
		//reset all constraints:


		resetconstraints(&(GetFoldDocument()->ct));
		
	}
	
}

void CFoldView::OnForceRestoreconstraints() 
{
	CFileDialog *filedialog;
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Constraint Files (*.con)|*.con|All Files|*.*||");

	
	
	if (filedialog->DoModal()==IDOK) {
		readconstraints(filedialog->GetPathName().GetBuffer(0),&(GetFoldDocument()->ct));
		

	}
	delete filedialog;
	
}

void CFoldView::OnForceSaveconstraints() 
{
	
	CFileDialog *filedialog;


	filedialog = new CFileDialog(FALSE,".con",NULL,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"Constraint Files (*.con)|*.con|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {


		outputconstraints(filedialog->GetPathName().GetBuffer(0),&(GetFoldDocument()->ct));

	

	}
	delete filedialog;

	
}

void CFoldView::OnForceSingestranded() 
{
	CSingleDialog *singledialog;

	singledialog = new CSingleDialog(&(GetFoldDocument()->ct),this);

	singledialog->Create(IDD_SINGLE_DIALOG,this);
	
}

/*This is unnecessary with the version 4 parameters
void CFoldView::OnAdvancedSuboptimalstructures() 
{
	//Open a Cadvanced dialog box to get info on suboptimal structure generation in
	//the dynamic programming algorithm
	CAdvanced *advancedialog;

	advancedialog = new CAdvanced(&dynnumber,&dynpercent,&dynwindow,this);

	advancedialog->Create(IDD_ADVANCED,this);
	
}*/

void CFoldView::OnForceFmn() 
{
	CSingleDialog *singledialog;

	singledialog = new CSingleDialog(&(GetFoldDocument()->ct),this,2);

	singledialog->Create(IDD_SINGLE_DIALOG,this);
	
}

void CFoldView::OnChangeParameter() 
{
	UpdateData(TRUE);
	
}

void CFoldView::OnForceChemical() 
{
	CSingleDialog *singledialog;

	singledialog = new CSingleDialog(&(GetFoldDocument()->ct),this,3);

	singledialog->Create(IDD_SINGLE_DIALOG,this);
	
}





void CFoldView::OnForceNmrconstraints()
{
	NMRDialog *nmrdialog;

	nmrdialog = new NMRDialog(&(GetFoldDocument()->ct),this);
	nmrdialog->Create(IDD_NMR,this);

}

void CFoldView::OnForceForbidbasepairs169()
{
	//Forbid Basepairs:
	
	CPairDialog *pairdialog;

	pairdialog = new CPairDialog(&(GetFoldDocument()->ct),this,true);

	pairdialog->Create(IDD_PAIR_DIALOG,this);

}

void CFoldView::OnTemperature()
{
	//Allow the user to specify a new temperature.
	CTemp_Dialog *temp;
	temp=new CTemp_Dialog(&(GetFoldDocument()->T));

	temp->DoModal();
	delete temp;
	


}

void CFoldView::OnMaximumloop()
{
	//User has chosen to change the maximum loop size.
	CInternalDialog *inter;
	inter=new CInternalDialog(&(GetFoldDocument()->maximuminternalloopsize));

	inter->DoModal();
	delete inter;

}

void CFoldView::OnForceReadshapereactvity()
{
	//Use has chosen to read SHAPE reactivity data:
	CReadSHAPE *readshape;
	readshape = new CReadSHAPE(&GetFoldDocument()->ct,this);
	readshape->DoModal();
	delete readshape;
}

void CFoldView::OnForceMaximumpairingdistance()
{
	//Use has chosen to set/change the maximum pairing distance:
	CDistanceDialog *dist;
	dist = new CDistanceDialog(&GetFoldDocument()->ct,this);
	dist->DoModal();
	delete dist;
}

void CFoldView::OnReadShapeLinear()
{
	//Use has chosen to read SHAPE reactivity data:
	CLinearSHAPEDialog *readshape;
	readshape = new CLinearSHAPEDialog(&GetFoldDocument()->ct,this);
	readshape->DoModal();
	delete readshape;
}

void CFoldView::OnSubfold()
{
	
	//User has chosen the SubFold menu option	
	CMenu* menu = GetFoldDocument()->menuframe->GetMenu( );	

	subfold = !subfold;

	if (subfold) menu->CheckMenuItem(ID_SUBFOLD,MF_CHECKED);
	else menu->CheckMenuItem(ID_SUBFOLD,MF_UNCHECKED); 


}

