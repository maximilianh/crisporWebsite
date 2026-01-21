// DistanceDialog.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "DistanceDialog.h"


// CDistanceDialog dialog

IMPLEMENT_DYNAMIC(CDistanceDialog, CDialog)
CDistanceDialog::CDistanceDialog(structure *CT,CWnd* pParent /*=NULL*/)
	: CDialog(CDistanceDialog::IDD, pParent)
	, m_distance(0)
	
{
	ct=CT;
	m_distance = ct->GetPairingDistanceLimit();
	

}

CDistanceDialog::~CDistanceDialog()
{
}

void CDistanceDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_DISTANCE, m_distance);
	DDV_MinMaxInt(pDX, m_distance, 5, 100000);
	
}

BOOL CDistanceDialog::OnInitDialog() {
	CButton *button;

	if (ct->DistanceLimited()) {
		button = (CButton*) GetDlgItem( IDC_RADIOYES );
		button->SetCheck(1);
	}
	else {
		button = (CButton*) GetDlgItem( IDC_RADIONO );
		button->SetCheck(1);
	}
	return CDialog::OnInitDialog();
}


BEGIN_MESSAGE_MAP(CDistanceDialog, CDialog)
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
END_MESSAGE_MAP()


// CDistanceDialog message handlers



void CDistanceDialog::OnBnClickedOk()
{
	CButton *button;
	UpdateData(TRUE);
	
	//get the button state
	button = (CButton*) GetDlgItem( IDC_RADIOYES );
	ct->SetPairingDistance(m_distance); //= (bool) button->GetCheck();

	//ct->maxdistance=;



	OnOK();
}
