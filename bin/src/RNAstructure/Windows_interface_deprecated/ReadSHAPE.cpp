// ReadSHAPE.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "ReadSHAPE.h"


// CReadSHAPE dialog

IMPLEMENT_DYNAMIC(CReadSHAPE, CDialog)
CReadSHAPE::CReadSHAPE(structure *CT, CWnd* pParent /*=NULL*/)
	: CDialog(CReadSHAPE::IDD, pParent)
	, InputFilename(_T(""))
	, ForceSingleStranded(2)
	, ChemicalModification(1)
{
	ct = CT;
}

CReadSHAPE::~CReadSHAPE()
{
}

void CReadSHAPE::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_InputFile, InputFilename);
	DDX_Text(pDX, IDC_SSThreshold, ForceSingleStranded);
	DDX_Text(pDX, IDC_ModThreshold, ChemicalModification);
}


BEGIN_MESSAGE_MAP(CReadSHAPE, CDialog)
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
	ON_BN_CLICKED(IDC_SHAPEFILEBUTTON, OnBnClickedShapefilebutton)
END_MESSAGE_MAP()


// CReadSHAPE message handlers

void CReadSHAPE::OnBnClickedOk()
{
	// Read and parse the SHAPE reactivity input file
	
	UpdateData(TRUE);

	
	if (InputFilename=="") {
		//Trap the case where the file was not previously specified.
		AfxMessageBox( "The SHAPE Date Input File was not specified.", 
				MB_OK|MB_ICONHAND);
		return;
	}

	ct->ReadSHAPE(InputFilename.GetBuffer(40),ForceSingleStranded,ChemicalModification);



	OnOK();
}

void CReadSHAPE::OnBnClickedShapefilebutton()
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
