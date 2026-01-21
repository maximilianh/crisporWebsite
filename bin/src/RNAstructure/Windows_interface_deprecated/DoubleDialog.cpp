// DoubleDialog.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "DoubleDialog.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CDoubleDialog dialog


CDoubleDialog::CDoubleDialog(structure *CT, CWnd* pParent /*=NULL*/)
	: CDialog(CDoubleDialog::IDD, pParent)
{
	//{{AFX_DATA_INIT(CDoubleDialog)
	m_double = 0;
	//}}AFX_DATA_INIT

	ct = CT;
}


void CDoubleDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDoubleDialog)
	DDX_Text(pDX, IDC_NUMBER_DOUBLE, m_double);
	DDV_MinMaxInt(pDX, m_double, 0, 32000);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CDoubleDialog, CDialog)
	//{{AFX_MSG_MAP(CDoubleDialog)
	ON_EN_CHANGE(IDC_NUMBER_DOUBLE, OnChangeNumberDouble)
	ON_BN_CLICKED(IDOK2, OnOk2)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDoubleDialog message handlers

void CDoubleDialog::OnChangeNumberDouble() 
{
	UpdateData(TRUE);
	if (ct->allocated) {
		if (m_double > ct->GetSequenceLength()) {
			
			AfxMessageBox( "Nucleotide is past end of sequence!", 
			MB_OK|MB_ICONHAND);
			m_double = 0;
			UpdateData(FALSE);
		}
	}
	
}

void CDoubleDialog::OnCancel() 
{
	
	
	CDialog::OnCancel();
}

void CDoubleDialog::OnOK() 
{
	
if (m_double!=0) {	
	//(ct->ndbl)++;

	//if// (ct->ndbl>maxforce) {
	//	AfxMessageBox( "Too many double-stranded nucleotides specified!",
	//		MB_OK|MB_ICONHAND);
	//	(ct->ndbl)--;
	//	CDialog::OnOK();
			

	//}

	ct->AddDouble(m_double);

	
}
CDialog::OnOK();

}

void CDoubleDialog::OnOk2() 
{
	if (m_double!=0) {	
	//(ct->ndbl)++;

	//if (ct->ndbl>maxforce) {
	//	AfxMessageBox( "Too many double-stranded nucleotides specified!",
	//		MB_OK|MB_ICONHAND);
	//	(ct->ndbl)--;
	//	CDialog::OnOK();
			

	//}

		ct->AddDouble( m_double);

	
}

	m_double = 0;
	UpdateData(FALSE);
	
}
