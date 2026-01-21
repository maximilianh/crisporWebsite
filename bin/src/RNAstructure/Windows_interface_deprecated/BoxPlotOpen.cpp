// BoxPlotOpen.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "BoxPlotOpen.h"


// CBoxPlotOpen dialog

IMPLEMENT_DYNAMIC(CBoxPlotOpen, CDialog)
CBoxPlotOpen::CBoxPlotOpen(CWnd* pParent /*=NULL*/)
	: CDialog(CBoxPlotOpen::IDD, pParent)
	, pfsave(_T(""))
{
}

CBoxPlotOpen::CBoxPlotOpen(CString *Filename, CWnd* pParent /*=NULL*/)
	: CDialog(CBoxPlotOpen::IDD, pParent)
	, pfsave(_T(""))
{

	filename = Filename;
}

CBoxPlotOpen::~CBoxPlotOpen()
{
}

void CBoxPlotOpen::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_PFSAVE, pfsave);
}


BEGIN_MESSAGE_MAP(CBoxPlotOpen, CDialog)
	ON_BN_CLICKED(IDC_OPENSAVE, OnBnClickedOpensave)
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
END_MESSAGE_MAP()


// CBoxPlotOpen message handlers

void CBoxPlotOpen::OnBnClickedOpensave()
{
	// Open a file open dialog to get the name of the save file
	CFileDialog *filedialog;
	

	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Partition Function Save Files (*.pfs)|*.dps||");

	
	//filedialog->m_ofn.lpstrInitialDir=GetFoldDocument()->startpath;
	if (filedialog->DoModal()==IDOK) {
		
		pfsave=(filedialog->GetPathName()).GetBuffer(30);
		
		UpdateData(FALSE);
	}
}

void CBoxPlotOpen::OnBnClickedOk()
{
	UpdateData(TRUE);
	if(pfsave=="") {
		MessageBox("Please specify a partition function save file name.");
		return;
	}
	
	*filename = pfsave;	

	OnOK();
}
