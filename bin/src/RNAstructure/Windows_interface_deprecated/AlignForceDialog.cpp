// AlignForceDialog.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "AlignForceDialog.h"


// CAlignForceDialog dialog

IMPLEMENT_DYNAMIC(CAlignForceDialog, CDialog)
CAlignForceDialog::CAlignForceDialog(short **Alignforce, int Length1, int Length2, CWnd* pParent /*=NULL*/)
	: CDialog(CAlignForceDialog::IDD, pParent)
	, m_ns1(0)
	, m_ns2(0)
{
	alignforce = Alignforce;
	length1=Length1;
	length2=Length2;
}

CAlignForceDialog::~CAlignForceDialog()
{
}

void CAlignForceDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_NS1, m_ns1);
	DDX_Text(pDX, IDC_NS2, m_ns2);
}


BEGIN_MESSAGE_MAP(CAlignForceDialog, CDialog)
	ON_BN_CLICKED(IDR_OKANDOPEN, OnBnClickedOkandopen)
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
END_MESSAGE_MAP()


// CAlignForceDialog message handlers

void CAlignForceDialog::OnBnClickedOkandopen()
{
	//Read the values and save them
	UpdateData(TRUE);
	if (m_ns1<0||m_ns2<0||m_ns1>length1||m_ns2>length2){
		AfxMessageBox( "One of the nucleotides is invalid.", 
			MB_OKCANCEL   |MB_ICONEXCLAMATION   );

	}
	else {
		alignforce[0][m_ns1]=m_ns2;
		alignforce[1][m_ns2]=m_ns1;
	}
}

void CAlignForceDialog::OnBnClickedOk()
{
	OnBnClickedOkandopen();
	OnOK();
}
