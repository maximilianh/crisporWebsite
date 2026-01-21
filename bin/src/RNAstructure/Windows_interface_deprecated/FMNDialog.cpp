// FMNDialog.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "FMNDialog.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CFMNDialog dialog


CFMNDialog::CFMNDialog(structure *CT, CWnd* pParent /*=NULL*/)
	: CDialog(CFMNDialog::IDD, pParent)
{
	//{{AFX_DATA_INIT(CFMNDialog)
	m_fmn = 0;
	//}}AFX_DATA_INIT
	ct = CT;
}


void CFMNDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CFMNDialog)
	DDX_Text(pDX, IDC_FMN, m_fmn);
	DDV_MinMaxInt(pDX, m_fmn, 0, 32000);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CFMNDialog, CDialog)
	//{{AFX_MSG_MAP(CFMNDialog)
	ON_EN_CHANGE(IDC_FMN, OnChangeFmn)
	//}}AFX_MSG_MAP
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CFMNDialog message handlers

void CFMNDialog::OnChangeFmn() 
{
	UpdateData(TRUE);
	if (ct->allocated) {
		if (m_fmn > ct->GetSequenceLength()) {
			
			AfxMessageBox( "Nucleotide is past end of sequence!", 
			MB_OK|MB_ICONHAND);
			m_fmn = 0;
			UpdateData(FALSE);
		}

		if (ct->numseq[m_fmn]!=4) {
			AfxMessageBox( "Nucleotide is not a U!", 
			MB_OK|MB_ICONHAND);
			m_fmn = 0;
			UpdateData(FALSE);


		}
	}
	
}

void CFMNDialog::OnBnClickedOk() 
{
	if (m_fmn!=0) {
		//ct->ngu++;
		//if (ct->ngu>maxforce) {
		//	AfxMessageBox( "Too many FMN cleavage points specified!",
		//		MB_OK|MB_ICONHAND);
		//	(ct->ngu)--;
		//	CDialog::OnOK();
			

		//}
		ct->AddGUPair(m_fmn);//gu[ct->ngu]=;
	}
	
	CDialog::OnOK();
}

