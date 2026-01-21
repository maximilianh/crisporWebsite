// DynRefoldView.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "DynRefoldView.h"
#include "dyndoc.h"
#include <direct.h>
#include "../src/dynalign.h"
#include "error.h"
#include <iostream>
using namespace std;

// DynRefoldView

IMPLEMENT_DYNCREATE(DynRefoldView, CFormView)

DynRefoldView::DynRefoldView()
	: CFormView(DynRefoldView::IDD)
	, percent(0)
	, maxtracebacks(0)
	, window(0)
	, alignwindow(0)
	, savefile(_T(""))
	, ct1file(_T(""))
	, ct2file(_T(""))
	, alignmentfile(_T(""))
{

	alignallocated = false;
}

DynRefoldView::~DynRefoldView()
{
	short i;
	if (alignallocated) {

		
		for (i=0;i<maxtracebacks;i++)  delete[] align[i];
		delete[] align;

	}
}

void DynRefoldView::OnInitialUpdate() 
{
	CFormView::OnInitialUpdate();
	
	ResizeParentToFit();
	
}

void DynRefoldView::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_PERCENT, percent);
	DDV_MinMaxInt(pDX, percent, 0, 100);
	DDX_Text(pDX, IDC_NUMBER, maxtracebacks);
	DDX_Text(pDX, IDC_WINDOW, window);
	DDV_MinMaxShort(pDX, window, 0, 9999);
	DDX_Text(pDX, IDC_ALIGNWINDOW, alignwindow);
	DDV_MinMaxShort(pDX, alignwindow, 0, 9999);
	DDX_Text(pDX, IDC_SEQUENCENAME, savefile);
	DDX_Text(pDX, IDC_CTNAME, ct1file);
	DDX_Text(pDX, IDC_CTNAME6, ct2file);
	DDX_Text(pDX, IDC_CTNAME7, alignmentfile);
}

BEGIN_MESSAGE_MAP(DynRefoldView, CFormView)
	ON_BN_CLICKED(IDC_SAVEBUTTON, OnBnClickedSavebutton)
	ON_BN_CLICKED(IDC_CTBUTTON, OnBnClickedCtbutton)
	ON_BN_CLICKED(IDC_CTBUTTON2, OnBnClickedCtbutton2)
	ON_BN_CLICKED(IDC_ALIGNMENT, OnBnClickedAlignment)
	ON_BN_CLICKED(IDC_START, OnBnClickedStart)
END_MESSAGE_MAP()


// DynRefoldView diagnostics

#ifdef _DEBUG
void DynRefoldView::AssertValid() const
{
	CFormView::AssertValid();
}

void DynRefoldView::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}
#endif //_DEBUG


// DynRefoldView message handlers

void DynRefoldView::OnBnClickedSavebutton()
{
			
	int j;
	int i;
	CFileDialog *filedialog;

		
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Save Files (*.dsv)|*.dsv||");

	
	filedialog->m_ofn.lpstrInitialDir=((CDynDoc*) GetDocument())->startpath;
	if (filedialog->DoModal()==IDOK) {
		//strcpy(m_sequencename.GetBuffer(10),(filedialog->GetPathName()).GetBuffer(0));
		savefile=(filedialog->GetPathName()).GetBuffer(30);
		
		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(((CDynDoc*) GetDocument())->startpath,_MAX_PATH);


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
		strncpy((((CDynDoc*) GetDocument())->startpath),path.GetBuffer(1),i);
		*(((CDynDoc*) GetDocument())->startpath + i) ='\0';

		//get the standard window sizes to fill in
		ifstream sav(savefile.GetBuffer(30),ios::binary);

		read(&sav,&j);//get through the flag for modification
		read(&sav,&i);//read the length of the first sequence
		sav.close();

		if (i>1200) {
   			window=20;
			alignwindow=3;
		}
		else if (i>800) {
   			window=15;
			alignwindow=3; 
		}
		else if (i>500) {
   			window=11;
           alignwindow=3; 
		}
		else if ((i)>300) {
   			window=7;
			alignwindow=2; 
		}
		else if ((i)>120) {
   			window=5;
			alignwindow=1;
            
		}
		else if ((i)>50) {
   			window=3;
			alignwindow=1;
            
		}
		else {
			window=2;
			alignwindow=0;
		}

		maxtracebacks = 20;
		percent = 20;


		ct1size = i;

	}
	
	delete filedialog;

	UpdateData(FALSE);
}

void DynRefoldView::OnBnClickedCtbutton()
{
	//The user is specifying the CT file name explicitly
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".ct","",OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"CT Files (*.ct)|*.ct|All Files|*.*||");
	filedialog->m_ofn.lpstrInitialDir=((CDynDoc*) GetDocument())->startpath;
	if (filedialog->DoModal()==IDOK) {

		ct1file=(filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);		
		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(((CDynDoc*) GetDocument())->startpath,_MAX_PATH);

		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(((CDynDoc*) GetDocument())->startpath,_MAX_PATH);


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
		strncpy(((CDynDoc*) GetDocument())->startpath,path.GetBuffer(1),i);
		*((((CDynDoc*) GetDocument())->startpath) + i) ='\0';

	}
	delete filedialog;
}

void DynRefoldView::OnBnClickedCtbutton2()
{
	//The user is specifying the CT file name explicitly
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".ct","",OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"CT Files (*.ct)|*.ct|All Files|*.*||");
	filedialog->m_ofn.lpstrInitialDir=((CDynDoc*) GetDocument())->startpath;
	if (filedialog->DoModal()==IDOK) {

		ct2file=(filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);		
		//_getcwd(((CDynDoc*) GetDocument())->startpath,_MAX_PATH);
		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(((CDynDoc*) GetDocument())->startpath,_MAX_PATH);


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
		strncpy(((CDynDoc*) GetDocument())->startpath,path.GetBuffer(1),i);
		*(((CDynDoc*) GetDocument())->startpath + i) ='\0';

	}
	delete filedialog;
}

void DynRefoldView::OnBnClickedAlignment()
{
	//The user is specifying the CT file name explicitly
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".ali","",OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"Alignment Files (*.ali)|*.ali|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		alignmentfile=(filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);	
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
		strncpy(((CDynDoc*) GetDocument())->startpath,path.GetBuffer(1),i);
		*(((CDynDoc*) GetDocument())->startpath + i) ='\0';

	}
	delete filedialog;
}

void DynRefoldView::OnBnClickedStart()
{
	short i;

	if(savefile==""||ct1file==""||ct2file=="") {
		MessageBox("Please specify names for the savefile and both ct output files.");
		return;
	}

	//OK - do the prediction!
	CDialog* finished;
	UpdateData(TRUE);

	align = new short *[maxtracebacks];//maximum number of tracebacks and next line and below at delete
	for (i=0;i<maxtracebacks;i++)  align[i] = new short [ct1size+1];
	alignallocated = true;
	

	//if refold Dynalign returns 14, report a traceback error to the user:
	if (refolddynalign(savefile.GetBuffer(10), &((CDynDoc*) GetDocument())->ct, &((CDynDoc*) GetDocument())->ct2, align, maxtracebacks, window, alignwindow, percent)==14) {
		
		
		CError *error=new CError();
		error->DoModal();
		delete error;

	}
	
	(((CDynDoc*) GetDocument())->ct).ctout(ct1file.GetBuffer(30));
	(((CDynDoc*) GetDocument())->ct2).ctout(ct2file.GetBuffer(30));
	alignout(align,alignmentfile.GetBuffer(10),&((CDynDoc*) GetDocument())->ct,&((CDynDoc*) GetDocument())->ct2);
	finished = new CDialog(IDD_FINISHED,this);

	if(finished->DoModal()==IDOK) {
		((CRNAstructureApp*) ((CDynDoc*) GetDocument())->pMainFrame)->Draw(ct1file.GetBuffer(10));
		((CRNAstructureApp*) ((CDynDoc*) GetDocument())->pMainFrame)->Draw(ct2file.GetBuffer(10));

	}
	

	delete finished;
	
	((CDynDoc*) GetDocument())->Frame->SendMessage(WM_CLOSE);

}
