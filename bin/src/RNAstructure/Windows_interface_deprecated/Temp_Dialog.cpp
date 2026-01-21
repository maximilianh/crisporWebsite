// Temp_Dialog.cpp : implementation file
// Dialog to query user for folding temperature

#include "stdafx.h"
#include "RNAstructure.h"
#include "Temp_Dialog.h"


// CTemp_Dialog dialog


IMPLEMENT_DYNAMIC(CTemp_Dialog, CDialog)
CTemp_Dialog::CTemp_Dialog(float *temp, CWnd* pParent /*=NULL*/)
	: CDialog(IDD_TEMP_DIALOG,  pParent)
{

	if (temp!=NULL) {
		m_temp = *temp;
		T = temp;
	}
}

CTemp_Dialog::~CTemp_Dialog()
{
}

void CTemp_Dialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_TEMPERATURE, m_temp);
}

BOOL CTemp_Dialog::OnInitDialog()
{
	CDialog::OnInitDialog();
	return TRUE;  // return TRUE  unless you set the focus to a control
}

BEGIN_MESSAGE_MAP(CTemp_Dialog, CDialog)
	ON_BN_CLICKED(IDOK, OnButtonOK)
	ON_BN_CLICKED(IDCANCEL, OnButtonCancel)
END_MESSAGE_MAP()





// CTemp_Dialog message handlers

void CTemp_Dialog::OnButtonOK()
{
	UpdateData(TRUE);
	*T = m_temp;
	CDialog::OnOK();
	
}

void CTemp_Dialog::OnButtonCancel()
{
	CDialog::OnCancel();
	
}
