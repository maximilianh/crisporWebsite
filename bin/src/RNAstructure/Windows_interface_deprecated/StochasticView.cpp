// C:\Documents and Settings\dhm\My Documents\RNAstructure\Windows_interface_deprecated\StochasticView.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "../src/stochastic.h"
#include "stochasticview.h"
#include <direct.h>

//This is the backend functions, set up to run in a seperate thread
UINT StochasticProc( LPVOID pParam ){    
	//For now, just do some things with compile-time definitions


	CStochasticObject* pObject = (CStochasticObject*) pParam;
	CStochasticView* pView = (CStochasticView*) pObject->parent; 
	
	structure ct;
	stochastic(&ct, pObject->savefile, pObject->m_structures,pObject->m_seed,pObject->progress);
	ct.ctout(pObject->ctoutfile);
		
      
	::PostMessage(pView->m_hWnd,ID_FOLDDONE,0,0);
	

	return 0;   // thread completed successfully
}


// CStochasticView dialog

IMPLEMENT_DYNCREATE(CStochasticView, CFormView)
CStochasticView::CStochasticView(CWnd* pParent /*=NULL*/)
	: CFormView(CStochasticView::IDD)
	, m_structures(1000)
	, m_seed(1234)
	, m_pfs(_T(""))
	, m_ct(_T(""))
{
	started = false;
}

CStochasticView::~CStochasticView()
{
}

void CStochasticView::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_STRUCTURES, m_structures);
	DDV_MinMaxInt(pDX, m_structures, 1, 32000000);
	DDX_Text(pDX, IDC_SEED, m_seed);
	DDV_MinMaxInt(pDX, m_seed, 1, 32000000);
	DDX_Text(pDX, IDC_SEQUENCENAME, m_pfs);
	DDX_Text(pDX, IDC_CTNAME, m_ct);
}


BEGIN_MESSAGE_MAP(CStochasticView, CFormView)
	ON_MESSAGE(ID_FOLDDONE, DoneFolding)
	ON_BN_CLICKED(IDC_PARTITION, OnBnClickedPartition)
	ON_BN_CLICKED(IDC_CTBUTTON, OnBnClickedCtbutton)
	ON_BN_CLICKED(IDC_START, OnBnClickedStart)
END_MESSAGE_MAP()


void CStochasticView::OnUpdate(CView*, LPARAM, CObject*)
{

	ResizeParentToFit(FALSE);
	ResizeParentToFit();
	

}

// CStochasticView message handlers

void CStochasticView::OnBnClickedPartition()
{
	//Open the standard dialog to pick a sequence
	char *ctname;
	short int i;
	CFileDialog *filedialog;
	
	

	
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Partition Function Save Files (*.pfs)|*.pfs||");

	
	filedialog->m_ofn.lpstrInitialDir=GetFoldDocument()->startpath;
	if (filedialog->DoModal()==IDOK) {
		//strcpy(m_sequencename.GetBuffer(10),(filedialog->GetPathName()).GetBuffer(0));
		m_pfs=(filedialog->GetPathName()).GetBuffer(30);
		i = m_pfs.GetLength();
		
		ctname = new char[i+4];//allocate enough space so that 
														//three characters can be added 
														//to the name if necessary
		strcpy(ctname,m_pfs.GetBuffer(10));
		//count the characters to the .
		
		while(i>=0){
			
			if (ctname[i]=='.') break;
			i--;
		}
		if (i==0) i = m_pfs.GetLength();
		strcpy(ctname+i+1,"ct\0");
		m_ct=ctname;
		
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

void CStochasticView::OnBnClickedCtbutton()
{
	//The user is specifying the CT file name explicitly
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".ct",m_ct,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"CT Files (*.ct)|*.ct|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		m_ct=(filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);		

	}
	delete filedialog;
}

void CStochasticView::OnBnClickedStart()
{

	//a trap to make sure the calculation is started more than once
	if (started) return;
	started=true;

	if(m_pfs=="") {
		MessageBox("Please specify a partion function save file name.");
		return;
	}
	UpdateData(TRUE);


	//TProgressDialog *progress;
	progress = new TProgressDialog();//parent,"Folding the RNA..."
	CRect *rect;
	rect = new CRect(10,40,210,120);
	
	if (GetFoldDocument()->ISRNA)
	progress->Create(NULL,"Folding the RNA...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);//(IDD_PROGRESSDIALOG);//
	else progress->Create(NULL,"Folding the DNA...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);
	delete rect;


	
	
	stochasticobject.parent=this;
	
	
	stochasticobject.m_structures = m_structures;
	stochasticobject.m_seed = m_seed;
	stochasticobject.ctoutfile = m_ct.GetBuffer(10);
	stochasticobject.savefile = m_pfs.GetBuffer(10);
	stochasticobject.progress = progress;
	

	
	
	AfxBeginThread(StochasticProc,&stochasticobject);
	

	

}

LRESULT CStochasticView::DoneFolding(WPARAM wParam, LPARAM lParam) {
	CDialog* finished;
	

	
	//offer to display the predicted structures if this is not subfolding
	finished = new CDialog(IDD_FINISHED,this);

	if(finished->DoModal()==IDOK) ((CRNAstructureApp*) GetFoldDocument()->pMainFrame)->Draw(m_ct.GetBuffer(10));
	

	delete finished;
	

	

	progress->SendMessage (WM_CLOSE);
	
	GetFoldDocument()->Frame->SendMessage(WM_CLOSE);
	
	return 0;


}


CFoldDoc *CStochasticView::GetFoldDocument() {
	
	return ((CFoldDoc*) GetDocument());	

}
