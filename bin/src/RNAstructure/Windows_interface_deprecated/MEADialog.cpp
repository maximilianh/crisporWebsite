// MEADialog.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "MEADialog.h"
#include "../src/MaxExpect.h"



//This is the backend functions, set up to run in a seperate thread
UINT MEAProc( LPVOID pParam ){    
	
	CMEAObject* pObject = (CMEAObject*) pParam;
	CMEADialog* pView = (CMEADialog*) pObject->parent; 
	
	structure ct;

	//Call the maximum expected accuracy algorithm
	bpMatch(&ct, pObject->savefile, pObject->gamma, pObject->percent, pObject->structures, pObject->window, pObject->progress);
	ct.ctout(pObject->ctoutfile);
		   
	::PostMessage(pView->m_hWnd,ID_FOLDDONE,0,0);
	

	return 0;   // thread completed successfully
}

// CMEADialog dialog

IMPLEMENT_DYNCREATE(CMEADialog, CFormView)

CMEADialog::CMEADialog(CWnd* pParent /*=NULL*/)
	: CFormView(CMEADialog::IDD)
	, m_SEQUENCE(_T(""))
	, m_CT(_T(""))
	, m_percent(50)
	, m_structures(1000)
	, m_window(5)
	, m_gamma(1)
{

	started = false;

}

CMEADialog::~CMEADialog()
{
}

void CMEADialog::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_SEQUENCENAME, m_SEQUENCE);
	DDX_Text(pDX, IDC_CTNAME, m_CT);
	DDX_Text(pDX, IDC_PERCENT, m_percent);
	DDX_Text(pDX, IDC_STRUCTURES, m_structures);
	//DDV_MinMaxInt(pDX, m_structures, 1, 32000000);
	DDX_Text(pDX, IDC_WINDOW, m_window);
	//DDV_MinMaxInt(pDX, m_window, 0, 32000);
	DDX_Text(pDX, IDC_GAMMA, m_gamma);
	//DDV_MinMaxDouble(pDX, m_gamma, 0, 32000);
}


BEGIN_MESSAGE_MAP(CMEADialog, CFormView)
	ON_MESSAGE(ID_FOLDDONE, DoneFolding)
	ON_BN_CLICKED(IDC_PARTITION, &CMEADialog::OnBnClickedPartition)
	ON_BN_CLICKED(IDC_CTBUTTON, &CMEADialog::OnBnClickedCtbutton)
	ON_BN_CLICKED(IDC_START, &CMEADialog::OnBnClickedStart)
END_MESSAGE_MAP()

// CPFFormView diagnostics

#ifdef _DEBUG
void CMEADialog::AssertValid() const
{
	CFormView::AssertValid();
}

void CMEADialog::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}
#endif //_DEBUG

void CMEADialog::OnUpdate(CView*, LPARAM, CObject*)
{

	ResizeParentToFit(FALSE);
	ResizeParentToFit();
	

}

// CMEADialog message handlers

void CMEADialog::OnBnClickedPartition()
{
	//Get the name of the partition function save file and post it to the view
		//Open the standard dialog to pick a sequence
	char *ctname;
	short int i;
	CFileDialog *filedialog;
	
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Partition Function Save Files (*.pfs)|*.pfs||");

	
	filedialog->m_ofn.lpstrInitialDir=GetFoldDocument()->startpath;
	if (filedialog->DoModal()==IDOK) {
		//strcpy(m_sequencename.GetBuffer(10),(filedialog->GetPathName()).GetBuffer(0));
		m_SEQUENCE=(filedialog->GetPathName()).GetBuffer(30);
		i = m_SEQUENCE.GetLength();
		
		ctname = new char[i+4];//allocate enough space so that 
														//three characters can be added 
														//to the name if necessary
		strcpy(ctname,m_SEQUENCE.GetBuffer(10));
		//count the characters to the .
		
		while(i>=0){
			
			if (ctname[i]=='.') break;
			i--;
		}
		if (i==0) i = m_SEQUENCE.GetLength();
		strcpy(ctname+i+1,"ct\0");
		m_CT=ctname;
		
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

void CMEADialog::OnBnClickedCtbutton()
{
	//The user is specifying the CT file name explicitly
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".ct",m_CT,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"CT Files (*.ct)|*.ct|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		m_CT=(filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);		

	}
	delete filedialog;
}

void CMEADialog::OnBnClickedStart()
{
	//OK, User clicked start, so do the calculation and then clean up
		
	//a trap to make sure the calculation is started more than once
	if (started) return;
	started=true;

	//Make sure the user specified a .pfs file:
	if(m_SEQUENCE=="") {
		MessageBox("Please specify a partion function save file name.");
		return;
	}
	UpdateData(TRUE);

	if(m_window<0) {
		MessageBox("Please specify a window that is greater than or equal to zero.");
		return;
	}
	if(m_structures<=0) {
		MessageBox("Please specify a number of structures that is greater than zero.");
		return;
	}
	if(m_percent<0) {
		MessageBox("Please enter a percent that is greater than or qual to 0.");
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

	
	meaobject.parent=this;
	
	meaobject.ctoutfile = m_CT.GetBuffer(10);
	meaobject.savefile = m_SEQUENCE.GetBuffer(10);
	meaobject.progress = progress;
	meaobject.window=m_window;
	meaobject.gamma=m_gamma;
	meaobject.structures=m_structures;
	meaobject.percent=m_percent;
	

	AfxBeginThread(MEAProc,&meaobject);
	

	
}


//This function is called by the backend thread
LRESULT CMEADialog::DoneFolding(WPARAM wParam, LPARAM lParam) {
	CDialog* finished;
	

	
	//offer to display the predicted structures if this is not subfolding
	finished = new CDialog(IDD_FINISHED,this);

	if(finished->DoModal()==IDOK) ((CRNAstructureApp*) GetFoldDocument()->pMainFrame)->Draw(m_CT.GetBuffer(10));
	

	delete finished;

	progress->SendMessage (WM_CLOSE);
	
	GetFoldDocument()->Frame->SendMessage(WM_CLOSE);
	
	return 0;


}


//Get the underlying
CFoldDoc *CMEADialog::GetFoldDocument() {
	
	return ((CFoldDoc*) GetDocument());	

}

