// AllFoldView.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "AllFoldView.h"

#include "doubledialog.h"
#include "pairdialog.h"
#include "singledialog.h"
#include "../src/outputconstraints.h"
#include "nmrdialog.h"
#include "../src/alltrace.h"
#include "foldview.h"
#include "regionalnmr.h"
#include "microarraydialog.h"
#include "temp_Dialog.h"


#include <direct.h>



//This is the backend functions, set up to run in a seperate thread
UINT AllFoldProc( LPVOID pParam ){    
	CAllFoldObject* pObject = (CAllFoldObject*) pParam; 
	AllFoldView* pView = (AllFoldView*) pObject->parent; 
	
	alltrace(pObject->ct,pObject->data,pObject->m_percent,(short) pObject->m_abs,pObject->progress,pObject->savefile);

	

	pObject->ct->ctout(pObject->ctoutfile);

	
	::PostMessage(pView->m_hWnd,ID_FOLDDONE,0,0);
	
	

	return 0;   // thread completed successfully
}



// AllFoldView dialog

IMPLEMENT_DYNCREATE(AllFoldView, CFormView)
AllFoldView::AllFoldView(CWnd* pParent /*=NULL*/)
	: CFormView(AllFoldView::IDD)
	, m_save(FALSE)
	, m_percent(0)
	, m_sequencename(_T(""))
	, m_ctname(_T(""))
	, m_abs(0)
{



}

AllFoldView::~AllFoldView()
{
}

void AllFoldView::OnInitialUpdate() {

	if (GetFoldDocument()->savefile!=NULL) {
		if (*GetFoldDocument()->savefile) {
			m_save = TRUE;
			UpdateData(FALSE);

		}
	}
	ResizeParentToFit();

	//set ct->stacking to true so that stacking can be tracked
	//DHM removed 5/10/11 to be consistent with other programs.
	//GetFoldDocument()->ct.stacking=true;


}
void AllFoldView::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	DDX_Check(pDX, IDC_SAVECHECK, m_save);
	DDX_Text(pDX, IDC_PERCENT, m_percent);
	DDV_MinMaxShort(pDX, m_percent, 0, 100);
	DDX_Text(pDX, IDC_SEQUENCENAME, m_sequencename);
	DDX_Text(pDX, IDC_CTNAME, m_ctname);
	DDX_Text(pDX, IDC_NUMBER, m_abs);
}


BEGIN_MESSAGE_MAP(AllFoldView, CFormView)
	ON_BN_CLICKED(IDC_SEQUENCEBUTTON, OnBnClickedSequencebutton)
	ON_BN_CLICKED(IDC_CTBUTTON, OnBnClickedCtbutton)
	
	ON_BN_CLICKED(IDC_START, OnBnClickedStart)
	ON_COMMAND(ID_FORCE_BASEPAIR, OnForceBasepair)
	ON_COMMAND(ID_FORCE_CHEMICALMODIFICATION, OnForceChemicalmodification)
	ON_COMMAND(ID_FORCE_DOUBLESTRANDED, OnForceDoublestranded)
	ON_COMMAND(ID_FORCE_FMN, OnForceFmn)
	ON_COMMAND(ID_FORCE_NMRCONSTRAINTS, OnForceNmrconstraints)
	ON_COMMAND(ID_FORCE_SINGESTRANDED, OnForceSingestranded)
	ON_COMMAND(ID_FORCE_FORBIDBASEPAIRS, OnForceForbidbasepairs169)
	ON_COMMAND(ID_FORCE_CURRENT, OnForceCurrent)
	ON_COMMAND(ID_FORCE_RESET, OnForceReset)
	ON_COMMAND(ID_FORCE_SAVECONSTRAINTS, OnForceSaveconstraints)
	ON_COMMAND(ID_FORCE_RESTORECONSTRAINTS, OnForceRestoreconstraints)
	ON_MESSAGE(ID_FOLDDONE, DoneFolding)
	ON_COMMAND(ID_FORCE_REGIONALNMRCONSTRAINTS, OnForceRegionalnmrconstraints)
	ON_COMMAND(ID_FORCE_MICROARRAYCONSTRAINTS, OnForceMicroarrayconstraints)
	ON_COMMAND(ID_TEMPERATURE, OnTemperature)
END_MESSAGE_MAP()


// AllFoldView message handlers

void AllFoldView::OnBnClickedSequencebutton()
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

CFoldDoc *AllFoldView::GetFoldDocument() {
	
	return ((CFoldDoc*) GetDocument());	

}

//call this function when a new sequence is chosen to set the default window size;
void AllFoldView::NewSequence() 
{

	
	
	//load the ct file:
	if((GetFoldDocument()->ct.openseq(m_sequencename.GetBuffer(10)))==0) {
   	//errmsg(5,5);
   }
   if ((GetFoldDocument()->ct.GetSequenceLength())>1200) {
   			m_percent=5;
			m_abs = .25;
   }
   else if ((GetFoldDocument()->ct.GetSequenceLength())>800) {
   			m_percent=8;
			m_abs = .5;
            
   }
   else if ((GetFoldDocument()->ct.GetSequenceLength())>500) {
   			m_percent=10;
			m_abs = .75;
            
   }
   else if ((GetFoldDocument()->ct.GetSequenceLength())>300) {
   			m_percent=15;
			m_abs = 1;
            
   }
   else if ((GetFoldDocument()->ct.GetSequenceLength())>120) {
   			m_percent=20;
			m_abs = 1.5;
            
   }
   else if ((GetFoldDocument()->ct.GetSequenceLength())>50) {
   			m_percent=25;
			m_abs = 3;
            
   }
   else {
	   m_percent=50;
		m_abs = 10;
   }



	UpdateData(FALSE);
}


void AllFoldView::OnBnClickedCtbutton()
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



void AllFoldView::OnBnClickedStart()
{
	int i;
	UpdateData(TRUE);
	if(m_sequencename=="") {
		MessageBox("Please specify a sequence name.");
		return;
	}

	//check to see if the temperature has been changed.
	if (GetFoldDocument()->T<310||GetFoldDocument()->T>311) {

		//change the temperature from 310.15
		if (GetFoldDocument()->newtemp()==0) {
			//if newtemp returned zero, pass a warning to the user
			AfxMessageBox( "An enthalpy data file could not be found!\nTemperature of prediction will revert back to 37 degrees C.\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
				MB_OK|MB_ICONHAND);

		}

	}

	
	progress = new TProgressDialog();//parent,"Folding the RNA..."
	CRect *rect;
	rect = new CRect(10,40,210,120);
	
	if (GetFoldDocument()->ISRNA)
	progress->Create(NULL,"Finding All RNA Tracebacks...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);//(IDD_PROGRESSDIALOG);//
	else progress->Create(NULL,"Finding All DNA Tracebacks...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);
	delete rect;


	foldobject.parent=this;
	foldobject.ct=&(GetFoldDocument()->ct);
	foldobject.data=&GetFoldDocument()->data;

	foldobject.m_percent = (int) m_percent;
	foldobject.m_abs = (int)(m_abs*10);
	foldobject.progress = progress;
	foldobject.ctoutfile = m_ctname.GetBuffer(10);
	

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

	

	
	AfxBeginThread(AllFoldProc,&foldobject);
}

void AllFoldView::OnForceBasepair()
{
	CPairDialog *pairdialog;

	pairdialog = new CPairDialog(&(GetFoldDocument()->ct),this);

	pairdialog->Create(IDD_PAIR_DIALOG,this);
}

void AllFoldView::OnForceChemicalmodification()
{
	CSingleDialog *singledialog;

	singledialog = new CSingleDialog(&(GetFoldDocument()->ct),this,3);

	singledialog->Create(IDD_SINGLE_DIALOG,this);
}

void AllFoldView::OnForceDoublestranded()
{
	CDoubleDialog *doubledialog;

	doubledialog = new CDoubleDialog(&(GetFoldDocument()->ct),this);

	doubledialog->Create(IDD_DOUBLE_DIALOG,this);
}

void AllFoldView::OnForceFmn()
{
	CSingleDialog *singledialog;

	singledialog = new CSingleDialog(&(GetFoldDocument()->ct),this,2);

	singledialog->Create(IDD_SINGLE_DIALOG,this);
}

void AllFoldView::OnForceNmrconstraints()
{
	NMRDialog *nmrdialog;

	nmrdialog = new NMRDialog(&(GetFoldDocument()->ct),this);
	nmrdialog->Create(IDD_NMR,this);
}

void AllFoldView::OnForceSingestranded()
{
	CSingleDialog *singledialog;

	singledialog = new CSingleDialog(&(GetFoldDocument()->ct),this);

	singledialog->Create(IDD_SINGLE_DIALOG,this);
}

void AllFoldView::OnForceForbidbasepairs169()
{
	//Forbid Basepairs:
	
	CPairDialog *pairdialog;

	pairdialog = new CPairDialog(&(GetFoldDocument()->ct),this,true);

	pairdialog->Create(IDD_PAIR_DIALOG,this);
}

void AllFoldView::OnForceCurrent()
{

	currentconstraints(&(GetFoldDocument()->ct));
}

void AllFoldView::OnForceReset()
{
	if (AfxMessageBox( "This will erase all constraints.\nContinue?", 
		MB_OKCANCEL   |MB_ICONQUESTION   )==IDOK) {
		//reset all constraints:
		resetconstraints(&(GetFoldDocument()->ct));
	}
}

void AllFoldView::OnForceSaveconstraints()
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

void AllFoldView::OnForceRestoreconstraints()
{
	CFileDialog *filedialog;
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Constraint Files (*.con)|*.con|All Files|*.*||");

	
	
	if (filedialog->DoModal()==IDOK) {
		readconstraints(filedialog->GetPathName().GetBuffer(0),&(GetFoldDocument()->ct));
		

	}
	delete filedialog;
}

LRESULT AllFoldView::DoneFolding(WPARAM wParam, LPARAM lParam) {
	CDialog* finished;
	//CMainFrame* frame = (CMainFrame*) Parent;

	if (m_save) delete foldobject.savefile;

	finished = new CDialog(IDD_FINISHED,this);

	if(finished->DoModal()==IDOK) ((CRNAstructureApp*) GetFoldDocument()->pMainFrame)->Draw(m_ctname.GetBuffer(10));
	

	delete finished;

	

	progress->SendMessage (WM_CLOSE);
	
	GetFoldDocument()->Frame->SendMessage(WM_CLOSE);
	
	return 0;


}

void AllFoldView::OnForceRegionalnmrconstraints()
{
	CRegionalNMR *nmrdialog;

	nmrdialog = new CRegionalNMR(&(GetFoldDocument()->ct),this);
	nmrdialog->Create(CRegionalNMR::IDD,this);
}

void AllFoldView::OnForceMicroarrayconstraints()
{
	CMicroarrayDialog *microdialog;
	microdialog = new CMicroarrayDialog(&(GetFoldDocument()->ct),this);
	microdialog->Create(CMicroarrayDialog::IDD,this);
}

void AllFoldView::OnTemperature()
{
	//Allow the user to specify a new temperature.
	CTemp_Dialog *temp;
	temp=new CTemp_Dialog(&(GetFoldDocument()->T));

	temp->DoModal();
	delete temp;
	


}
