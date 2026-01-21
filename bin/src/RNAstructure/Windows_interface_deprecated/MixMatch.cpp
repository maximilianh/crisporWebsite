// MixMatch.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "MixMatch.h"
#include "../src/mm.cpp"

#include <direct.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

UINT MixMatchProc( LPVOID pParam ){    
	CMixMatchObject* pObject = (CMixMatchObject*) pParam;
	CMixMatch* pView = (CMixMatch*) pObject->parent; 
    int i;
	structure ct1,ct2;
	
	openct(&ct1,pObject->ctinfile);

   

	ct2.allocate(ct1.numofbases);

	//copy over the historical numbering and characters to ct2
	for (i=1;i<=ct1.numofbases;i++) {
		ct2.nucs[i]=ct1.nucs[i];
		ct2.hnumber[i] = ct1.hnumber[i];

	}


	mixmatch(&ct1,&ct2,pObject->data,pObject->modfile, pObject->progress);

 

	ctout(&ct2,pObject->ctoutfile);

      

	
	
	::PostMessage(pView->m_hWnd,ID_MIXMATCHDONE,0,0);
	

	return 0;   // thread completed successfully
}


/////////////////////////////////////////////////////////////////////////////
// CMixMatch dialog


CMixMatch::CMixMatch()
	: CFormView(CMixMatch::IDD)
{
	//{{AFX_DATA_INIT(CMixMatch)
	m_ctinname = _T("");
	m_ctoutname = _T("");
	m_modname = _T("");
	//}}AFX_DATA_INIT
	
}

IMPLEMENT_DYNCREATE(CMixMatch, CFormView)
void CMixMatch::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CMixMatch)
	DDX_Text(pDX, IDC_CTINNAME, m_ctinname);
	DDX_Text(pDX, IDC_OUTCTNAME, m_ctoutname);
	DDX_Text(pDX, IDC_MODNAME, m_modname);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CMixMatch, CFormView)
	ON_MESSAGE(ID_MIXMATCHDONE, Done)
	//{{AFX_MSG_MAP(CMixMatch)
	ON_BN_CLICKED(IDC_INCT, OnInct)
	ON_BN_CLICKED(IDC_MOD, OnMod)
	ON_BN_CLICKED(IDC_CTOUT, OnCtout)
	ON_BN_CLICKED(IDOK, OnOK)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CMixMatch message handlers

void CMixMatch::OnInitialUpdate() {

	ResizeParentToFit();
}

void CMixMatch::OnInct() 
{
	CFileDialog *filedialog;
	
	

	
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"CT Files (*.ct)|*.ct||");

	
	filedialog->m_ofn.lpstrInitialDir=GetFoldDocument()->startpath;
	if (filedialog->DoModal()==IDOK) {
		//strcpy(m_sequencename.GetBuffer(10),(filedialog->GetPathName()).GetBuffer(0));
		m_ctinname=(filedialog->GetPathName()).GetBuffer(30);
		
		UpdateData(FALSE);
		
		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(GetFoldDocument()->startpath,_MAX_PATH);

		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(GetFoldDocument()->startpath,_MAX_PATH);
		CString path;
		path = filedialog->GetPathName();
		int i = path.GetLength();
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

void CMixMatch::OnMod() 
{
	CFileDialog *filedialog;
	
	

	
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"mod Files (*.mod)|*.mod||");

	
	filedialog->m_ofn.lpstrInitialDir=GetFoldDocument()->startpath;
	if (filedialog->DoModal()==IDOK) {
		//strcpy(m_sequencename.GetBuffer(10),(filedialog->GetPathName()).GetBuffer(0));
		m_modname=(filedialog->GetPathName()).GetBuffer(30);
		
		UpdateData(FALSE);
		
		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(GetFoldDocument()->startpath,_MAX_PATH);

		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(GetFoldDocument()->startpath,_MAX_PATH);
		CString path;
		path = filedialog->GetPathName();
		int i = path.GetLength();
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

void CMixMatch::OnCtout() 
{
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".ct",NULL,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"CT Files (*.ct)|*.ct|All Files|*.*||");
	filedialog->m_ofn.lpstrInitialDir=GetFoldDocument()->startpath;
	if (filedialog->DoModal()==IDOK) {

		m_ctoutname=(filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);
		
		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(GetFoldDocument()->startpath,_MAX_PATH);

		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(GetFoldDocument()->startpath,_MAX_PATH);
		CString path;
		path = filedialog->GetPathName();
		int i = path.GetLength();
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

void CMixMatch::OnOK() 
{

	
	if (m_ctinname==""||m_ctoutname==""||m_modname=="") {

		MessageBox("Please specify a name for all three fields.");
		return;

	}

	//Set up a thread to do the mix and match
	progress = new TProgressDialog();//parent,"Folding the RNA..."
	CRect *rect;
	rect = new CRect(10,40,210,120);
	progress->Create(NULL,"Mix&Match in progress...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);
	
	object.progress = progress;
	object.ctinfile = m_ctinname.GetBuffer(0);
	object.data = &(GetFoldDocument()->data);
	object.parent = this;
	object.ctoutfile = m_ctoutname.GetBuffer(0);
	object.modfile = m_modname.GetBuffer(0);

	AfxBeginThread(MixMatchProc,&object);
	
}

CFoldDoc *CMixMatch::GetFoldDocument() {
	
	return ((CFoldDoc*) GetDocument());	

}

LRESULT CMixMatch::Done(WPARAM wParam, LPARAM lParam) {
	progress->SendMessage (WM_CLOSE);
	
	CDialog* finished;
	

	finished = new CDialog(IDD_FINISHED,this);

	if(finished->DoModal()==IDOK) ((CRNAstructureApp*) GetFoldDocument()->pMainFrame)->Draw(m_ctoutname.GetBuffer(10));
	

	delete finished;

	


	GetFoldDocument()->Frame->SendMessage(WM_CLOSE);

	return 0;
}