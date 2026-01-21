// LinearSHAPEDialog.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "LinearSHAPEDialog.h"


// CLinearSHAPEDialog dialog

IMPLEMENT_DYNAMIC(CLinearSHAPEDialog, CDialog)
CLinearSHAPEDialog::CLinearSHAPEDialog(structure *CT, CWnd* pParent /*=NULL*/)
	: CDialog(CLinearSHAPEDialog::IDD, pParent)
	, Slope(1.8)
	, InputFilename(_T(""))
	, intercept(-0.6)
{
	ct = CT;
}

CLinearSHAPEDialog::~CLinearSHAPEDialog()
{
}

void CLinearSHAPEDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_InputFile, InputFilename);
	DDX_Text(pDX, IDC_SLOPE, Slope);
	DDX_Text(pDX, IDC_INTERCEPT, intercept);
}


BEGIN_MESSAGE_MAP(CLinearSHAPEDialog, CDialog)
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
	ON_BN_CLICKED(IDCANCEL, OnBnClickedCancel)
	ON_BN_CLICKED(IDC_SHAPEFILEBUTTON, OnBnClickedShapefilebutton)
END_MESSAGE_MAP()


// CLinearSHAPEDialog message handlers

void CLinearSHAPEDialog::OnBnClickedOk()
{
	/// Read and parse the SHAPE reactivity input file
	
	UpdateData(TRUE);

	
	if (InputFilename=="") {
		//Trap the case where the file was not previously specified.
		AfxMessageBox( "The SHAPE Date Input File was not specified.", 
				MB_OK|MB_ICONHAND);
		return;
	}

	ct->SHAPEslope=Slope*conversionfactor;
	ct->SHAPEintercept=intercept*conversionfactor;
	ct->ReadSHAPE(InputFilename.GetBuffer(40));
	


	OnOK();
}

void CLinearSHAPEDialog::OnBnClickedCancel()
{
	// TODO: Add your control notification handler code here
	OnCancel();
}

void CLinearSHAPEDialog::OnBnClickedShapefilebutton()
{
	//Get the input filename
	CFileDialog *filedialog;

	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"SHAPE Files (*.shape)|*.shape|All Files (*.*)|*.*||");

	
	if (filedialog->DoModal()==IDOK) {
		//strcpy(m_sequencename.GetBuffer(10),(filedialog->GetPathName()).GetBuffer(0));
		InputFilename=(filedialog->GetPathName()).GetBuffer(30);
	}
	UpdateData(FALSE);

	delete filedialog;
}
