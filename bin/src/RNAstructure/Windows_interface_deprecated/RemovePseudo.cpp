// RemovePseudo.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "RemovePseudo.h"
#include "../src/alltrace.h"
#include "temp_Dialog.h"


#include <direct.h>

//This is the backend functions, set up to run in a seperate thread
UINT BPProc( LPVOID pParam ){    
	CAllFoldObject* pObject = (CAllFoldObject*) pParam; 
	CRemovePseudo* pView = (CRemovePseudo*) pObject->parent; 
	
	
	int i,j,structures,tracebackstatus;
	structure *tempct;


	//allocate tempct:
	tempct = new structure(2);
	tempct->allocate(pObject->ct->GetSequenceLength());
	tempct->SetSequenceLabel("temp\n");
	//strcpy(tempct->ctlabel[1],;
	//tempct->numofbases = pObject->ct->numofbases;

	for (i=1;i<=pObject->ct->GetSequenceLength();i++) {
		tempct->numseq[i]=pObject->ct->numseq[i];

	}
	//allocate a template of allowed pairs
	tempct->allocatetem();

	//Break pairs for each structure
	for (structures=1;structures<=pObject->ct->GetNumberofStructures();structures++) {

		//initialize all base pairs as unallowed:
		for (i=0;i<=pObject->ct->GetSequenceLength();i++) {
			for (j=i+1;j<=pObject->ct->GetSequenceLength();j++) {
    			tempct->tem[j][i] = false;
			}
	   }

		//now allow the pairs that are in the loaded ct:
		for (i=1;i<=pObject->ct->GetSequenceLength();i++) {
			if (pObject->ct->GetPair(i,structures)>i) {
				tempct->tem[pObject->ct->GetPair(i,structures)][i] = true;
			}
		}
			
		//tempct->numofstructures = 0; //strip the structure from the ct at this point
		while (tempct->GetNumberofStructures()>0) tempct->RemoveLastStructure();
		
		//Predict the secondary structures.
		tracebackstatus = 0;//no error tracking yet in alltrace
		alltrace(tempct,pObject->data,0,0,pObject->progress,NULL,false);
		//dynamic(tempct, pObject->data, 1, 0, 0, pObject->progress, false, NULL, 30);

		//copy the pairs back to ct
		for (i=1;i<=pObject->ct->GetSequenceLength();i++) {
			//if (tempct->GetPair(i,1)>i) {
				pObject->ct->SetPair(i,tempct->GetPair(i,1),structures);
			//}

		}

		//Also copy the energy back to ct
		pObject->ct->SetEnergy(structures,tempct->GetEnergy(1));
		
		

	}
	delete tempct;
	
	
	//alltrace(pObject->ct,pObject->data,0,0,pObject->progress,NULL);

	
	pObject->ct->ctout (pObject->ctoutfile);

	

	
	::PostMessage(pView->m_hWnd,ID_FOLDDONE,0,0);
	
	

	return 0;   // thread completed successfully
}
// CRemovePseudo

IMPLEMENT_DYNCREATE(CRemovePseudo, CFormView)

void CRemovePseudo::OnInitialUpdate() {

	
	ResizeParentToFit();

	


}

CRemovePseudo::CRemovePseudo()
	: CFormView(CRemovePseudo::IDD)
{
}

CRemovePseudo::~CRemovePseudo()
{
}

void CRemovePseudo::DoDataExchange(CDataExchange* pDX)
{
	
	CFormView::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_CTINNAME, m_ctinname);
	DDX_Text(pDX, IDC_CTOUTNAME, m_ctoutname);
}

BEGIN_MESSAGE_MAP(CRemovePseudo, CFormView)
	ON_BN_CLICKED(IDC_CTIN, OnBnClickedCtin)
	ON_BN_CLICKED(IDC_CTBUTTON, OnBnClickedCtbutton)
	ON_BN_CLICKED(IDC_START, OnBnClickedStart)
	ON_MESSAGE(ID_FOLDDONE, DoneFolding)
	ON_COMMAND(ID_TEMPERATURE, OnTemperature)

END_MESSAGE_MAP()


// CRemovePseudo diagnostics

#ifdef _DEBUG
void CRemovePseudo::AssertValid() const
{
	CFormView::AssertValid();
}

void CRemovePseudo::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}
#endif //_DEBUG


CFoldDoc *CRemovePseudo::GetFoldDocument() {
	
	return ((CFoldDoc*) GetDocument());	

}


// CRemovePseudo message handlers

void CRemovePseudo::OnBnClickedCtin()
{

	//user is inputing sequence name
	char *ctname;
	short int i;
	CFileDialog *filedialog;
	
	

	
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"CT Files (*.ct)|*.ct||");

	
	filedialog->m_ofn.lpstrInitialDir=GetFoldDocument()->startpath;
	if (filedialog->DoModal()==IDOK) {
		//strcpy(m_sequencename.GetBuffer(10),(filedialog->GetPathName()).GetBuffer(0));
		m_ctinname=(filedialog->GetPathName()).GetBuffer(30);
		i = m_ctinname.GetLength();
		
		ctname = new char[i+15];//allocate enough space for the seq name plus the extra characters
		strcpy(ctname,m_ctinname.GetBuffer(10));
		//count the characters to the .
		
		while(i>=0){
			
			if (ctname[i]=='.') break;
			i--;
			
		}
		i--;
		if (i==0) i = m_ctinname.GetLength();
		strcpy(ctname+i+1,"_no_psuedo.ct\0");
		m_ctoutname=ctname;
		
		delete[] ctname;//fix this?
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

		

	}
	delete filedialog;
}

void CRemovePseudo::OnBnClickedCtbutton()
{
	//The user is specifying the output CT file name explicitly
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".ct",m_ctoutname,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"CT Files (*.ct)|*.ct|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		m_ctoutname=(filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);		

	}
	delete filedialog;
}

void CRemovePseudo::OnBnClickedStart()
{
	short i,j;
	UpdateData(TRUE);
	if(m_ctinname=="") {
		MessageBox("Please specify a input ct name.");
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

	GetFoldDocument()->ct.openct(m_ctinname.GetBuffer(0));
	(GetFoldDocument()->ct).allocatetem();

	//add a new line to the end of the title
	for (i=1;i<=(GetFoldDocument()->ct).GetNumberofStructures();i++)
		//strcat((GetFoldDocument()->ct).ctlabel[i],"\n");
		
		(GetFoldDocument()->ct).SetSequenceLabel("\n");
	//Move the following lines to the seperate thread and add multiple structure support, 9/2/09:
	//initialize all base pairs as unallowed:
	//for (i=0;i<=(GetFoldDocument()->ct).numofbases;i++) {
	//	for (j=i+1;j<=(GetFoldDocument()->ct).numofbases;j++) {
    //		(GetFoldDocument()->ct).tem[j][i] = false;
	//	}
   //}

	//now allow the pairs that are in the loaded ct:
	//for (i=0;i<=(GetFoldDocument()->ct).numofbases;i++) {
	//	if ((GetFoldDocument()->ct).basepr[1][i]>i) {
	//		(GetFoldDocument()->ct).tem[(GetFoldDocument()->ct).basepr[1][i]][i] = true;
	//	}
	//}
		
	//(GetFoldDocument()->ct).numofstructures = 0; //strip the structure from the ct at this point
	foldobject.parent=this;
	foldobject.ct=&(GetFoldDocument()->ct);
	foldobject.data=&GetFoldDocument()->data;

	
	foldobject.progress = progress;
	foldobject.ctoutfile = m_ctoutname.GetBuffer(10);
	

	
	foldobject.savefile = 0;
		
	
	

	

	
	AfxBeginThread(BPProc,&foldobject);
}

void CRemovePseudo::OnTemperature()
{
	//Allow the user to specify a new temperature.
	CTemp_Dialog *temp;
	temp=new CTemp_Dialog(&(GetFoldDocument()->T));

	temp->DoModal();
	delete temp;
	


}

LRESULT CRemovePseudo::DoneFolding(WPARAM wParam, LPARAM lParam) {
	CDialog* finished;
	//CMainFrame* frame = (CMainFrame*) Parent;

	

	finished = new CDialog(IDD_FINISHED,this);

	if(finished->DoModal()==IDOK) ((CRNAstructureApp*) GetFoldDocument()->pMainFrame)->Draw(m_ctoutname.GetBuffer(10));
	

	delete finished;

	

	progress->SendMessage (WM_CLOSE);
	
	GetFoldDocument()->Frame->SendMessage(WM_CLOSE);
	
	return 0;


}
