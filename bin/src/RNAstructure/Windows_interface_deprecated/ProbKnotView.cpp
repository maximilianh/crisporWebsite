// ProbKnotView.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "ProbKnotView.h"

#include "..\src\pfunction.h"
#include "..\src\structure.h"
#include "..\src\probknot.h"

using namespace std;

//This is the backend functions, set up to run in a seperate thread
UINT PKProc( LPVOID pParam ){    
	
	CPKObject* pObject = (CPKObject*) pParam;
	CProbKnotView* pView = (CProbKnotView*) pObject->parent; 

	pfunctionclass *w,*v,*wmb,*wmbl,*wcoax,*wl;
	forceclass* fce;
	PFPRECISION *w5, *w3;
	bool *lfce, *mod;
	short vers;
	structure ct;

	pfdatatable *data;

	//allocate the ct file by reading the save file:
	ifstream sav(pObject->savefile,ios::binary);

	
	read(&sav,&(vers));//read the version of the save file
		//right now there is no infrastructure to indicate the wrong version is being read. 
		//This should be changed in the future...

	int length;
	read(&sav,&(length));
	
	sav.close();
	//allocate everything
	
	ct.allocate(length);

	w = new pfunctionclass(ct.GetSequenceLength());
	v = new pfunctionclass(ct.GetSequenceLength());
	wmb = new pfunctionclass(ct.GetSequenceLength());
	wmbl = new pfunctionclass(ct.GetSequenceLength());
	wcoax = new pfunctionclass(ct.GetSequenceLength());
	wl = new pfunctionclass(ct.GetSequenceLength());
	wlc = new pfunctionclass(ct.GetSequenceLength());
	fce = new forceclass(ct.GetSequenceLength());

	w5 = new PFPRECISION [ct.GetSequenceLength()+1];
	w3 = new PFPRECISION [ct.GetSequenceLength()+2];

	lfce = new bool [2*ct.GetSequenceLength()+1];
    mod = new bool [2*ct.GetSequenceLength()+1];

	data = new pfdatatable();

	

	//load all the data from the pfsavefile:
	readpfsave(pObject->savefile, &ct, w5, w3,v, w, wmb,wl, wlc, wmbl, wcoax, fce,&data->scaling,mod,lfce,data);


	//Call the maximum expected accuracy algorithm
	
	ProbKnotAssemble(v, w5, &ct, data, lfce, mod, data->scaling, fce, pObject->iterations, pObject->minlength );

	ct.ctout(pObject->ctoutfile);
		   
	::PostMessage(pView->m_hWnd,ID_FOLDDONE,0,0);

	delete w;
	delete v ;
	delete wmb ;
	delete wmbl ;
	delete wcoax ;
	delete wl ;
	delete fce ;

	delete[] w5 ;
	delete[] w3 ;

	delete[] lfce ;
    delete[] mod ;

	delete data ;
	

	return 0;   // thread completed successfully
}


// CProbKnotView

IMPLEMENT_DYNCREATE(CProbKnotView, CFormView)

CProbKnotView::CProbKnotView()
	: CFormView(CProbKnotView::IDD)
	, m_savefile(_T(""))
	, m_ctfilename(_T(""))
	, m_iterations(1)
	, m_ml(3)
{

	started = false;

}

CProbKnotView::~CProbKnotView()
{
}

void CProbKnotView::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_SAVENAME, m_savefile);
	DDX_Text(pDX, IDC_CTNAME, m_ctfilename);
	DDX_Text(pDX, IDC_ITERATIONS, m_iterations);
	DDV_MinMaxInt(pDX, m_iterations, 1, 32000);
	DDX_Text(pDX, IDC_ML, m_ml);
	DDV_MinMaxInt(pDX, m_ml, 1, 32000);
}

BEGIN_MESSAGE_MAP(CProbKnotView, CFormView)
	ON_MESSAGE(ID_FOLDDONE, DoneFolding)
	ON_BN_CLICKED(IDC_OPENSAVEBUTTON, &CProbKnotView::OnBnClickedOpensavebutton)
	ON_BN_CLICKED(IDC_CTBUTTON, &CProbKnotView::OnBnClickedCtbutton)
	ON_BN_CLICKED(IDC_START, &CProbKnotView::OnBnClickedStart)
END_MESSAGE_MAP()


// CProbKnotView diagnostics

#ifdef _DEBUG
void CProbKnotView::AssertValid() const
{
	CFormView::AssertValid();
}

#ifndef _WIN32_WCE
void CProbKnotView::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}
#endif
#endif //_DEBUG


// CProbKnotView message handlers


//User clicked the savefilename button.
void CProbKnotView::OnBnClickedOpensavebutton()
{
	//Open the standard dialog to pick a sequence
	char *tsavename;
	short int i;
	CFileDialog *filedialog;
	
	

	
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Partition Function Save Files (*.pfs)|*.pfs||");

	
	filedialog->m_ofn.lpstrInitialDir=GetFoldDocument()->startpath;
	if (filedialog->DoModal()==IDOK) {
		//strcpy(m_sequencename.GetBuffer(10),(filedialog->GetPathName()).GetBuffer(0));

		//open the save file to check file version:
		short vers;
		ifstream sav((filedialog->GetPathName()).GetBuffer(30),ios::binary);	
		read(&sav,&(vers));//read the version of the save file
		sav.close();
		if (vers!=pfsaveversion) {
			//the file version doesn't match the current version
			AfxMessageBox( "Error: This save file was created with a different version of RNAstructure.  \nPlease select another file.", 
			MB_OK|MB_ICONEXCLAMATION);
			delete filedialog;
			return;


		}

		m_savefile=(filedialog->GetPathName()).GetBuffer(30);
		i = m_savefile.GetLength();
		
		
		
		
		tsavename = new char[i+2];//allocate enough space so that 
														//two characters can be added 
														//to the name if necessary
		
		
		
		strcpy(tsavename,m_savefile.GetBuffer(10));
		//count the characters to the .
		
		while(i>=0){
			
			if (tsavename[i]=='.') break;
			i--;
		}
		if (i==0) i = m_savefile.GetLength();
		strcpy(tsavename+i+1,"ct\0");
		m_ctfilename=tsavename;
		
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

			
		

	}
	delete filedialog;


	

	
}


//User Clicked the Button to Change the Output Filename.
void CProbKnotView::OnBnClickedCtbutton()
{
	//The user is specifying the CT file name explicitly
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".ct",m_ctfilename,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"CT Files (*.ct)|*.ct|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		m_ctfilename=(filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);		

	}
	delete filedialog;
}

void CProbKnotView::OnBnClickedStart()
{
	
	//OK, User clicked start, so do the calculation and then clean up
		
	//a trap to make sure the calculation is started more than once
	if (started) return;
	started=true;

	//Make sure the user specified a .pfs file:
	if(m_savefile=="") {
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

	
	pkobject.parent=this;
	
	pkobject.ctoutfile = m_ctfilename.GetBuffer(10);
	pkobject.savefile = m_savefile.GetBuffer(10);
	pkobject.progress = progress;
	pkobject.minlength=m_ml;
	pkobject.iterations=m_iterations;
	
	

	AfxBeginThread(PKProc,&pkobject);
}


//Get the underlying document
CFoldDoc *CProbKnotView::GetFoldDocument() {
	
	return ((CFoldDoc*) GetDocument());	

}


void CProbKnotView::OnUpdate(CView*, LPARAM, CObject*)
{

	ResizeParentToFit(FALSE);
	ResizeParentToFit();

	
	

}

//This function is called by the backend thread
LRESULT CProbKnotView::DoneFolding(WPARAM wParam, LPARAM lParam) {
	CDialog* finished;
	

	
	//offer to display the predicted structures if this is not subfolding
	finished = new CDialog(IDD_FINISHED,this);

	if(finished->DoModal()==IDOK) ((CRNAstructureApp*) GetFoldDocument()->pMainFrame)->Draw(m_ctfilename.GetBuffer(10));
	

	delete finished;

	progress->SendMessage (WM_CLOSE);
	
	GetFoldDocument()->Frame->SendMessage(WM_CLOSE);
	
	return 0;


}