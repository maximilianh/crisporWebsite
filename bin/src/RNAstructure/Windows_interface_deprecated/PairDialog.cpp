// PairDialog.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "PairDialog.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CPairDialog dialog


CPairDialog::CPairDialog(structure *CT,CWnd* pParent,bool Isforbid /*=NULL*/)
	: CDialog(CPairDialog::IDD, pParent)
{
	//{{AFX_DATA_INIT(CPairDialog)
	m_base1 = 0;
	m_base2 = 0;
	m_length = 0;
	//}}AFX_DATA_INIT
	ct = CT;
	isforbid = Isforbid;
}


void CPairDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CPairDialog)
	DDX_Text(pDX, IDC_BASE1, m_base1);
	DDV_MinMaxInt(pDX, m_base1, 0, 32000);
	DDX_Text(pDX, IDC_BASE2, m_base2);
	DDV_MinMaxInt(pDX, m_base2, 0, 32000);
	DDX_Text(pDX, IDC_LENGTH, m_length);
	DDV_MinMaxInt(pDX, m_length, 0, 32000);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CPairDialog, CDialog)
	//{{AFX_MSG_MAP(CPairDialog)
	ON_EN_CHANGE(IDC_BASE1, OnChangeBase1)
	ON_EN_CHANGE(IDC_BASE2, OnChangeBase2)
	ON_EN_CHANGE(IDC_LENGTH, OnChangeLength)
	ON_BN_CLICKED(IDOK2, OnOk2)
	//}}AFX_MSG_MAP
	//ON_BN_CLICKED(IDOK3, OnBnClickedOk3)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CPairDialog message handlers

void CPairDialog::OnChangeBase1() 
{
	
	
}

void CPairDialog::OnChangeBase2() 
{
	//UpdateData(TRUE);
	
	
}

void CPairDialog::OnChangeLength() 
{
	//UpdateData(TRUE);
	
	
}

void CPairDialog::OnOK() 
{

	getinfo();
	
	CDialog::OnOK();
}

void CPairDialog::OnOk2() 
{
	getinfo();
	

	m_base1=0;
	m_base2=0;
	m_length=0;
	UpdateData(FALSE);	
	
	
	
}

void CPairDialog::OnBnClickedOk3()
{
	// TODO: Add your control notification handler code here
}

void CPairDialog::getinfo() {

	short i;
	UpdateData(TRUE);

	if (m_base2<m_base1) {
		AfxMessageBox( "Nucleotide two should be 3' to nucleotide 1.", 
			MB_OK|MB_ICONHAND);
			m_base2 = 0;
			UpdateData(FALSE);
			return;
	}
	if (m_length>(m_base2-m_base1)/2) {
		AfxMessageBox( "Helix length is too long.", 
			MB_OK|MB_ICONHAND);
			m_length = 0;
			UpdateData(FALSE);
			return;
	}
	for (i=m_base1;i<m_base1+m_length;i++) {
	
		
		//if (!isforbid&&ct->npair==maxforce) {
		//	AfxMessageBox( "Too many paired nucleotides specified!",
		//		MB_OK|MB_ICONHAND);
		//	CDialog::OnOK();
		//}
		//else if (isforbid&&ct->nforbid==maxforce) {
		//	AfxMessageBox( "Too many forbidden pairs specified!",
		//		MB_OK|MB_ICONHAND);
		//	CDialog::OnOK();
		//
		//}
		if (isforbid) {
			//ct->forbid[ct->nforbid][0]=i;
			//ct->forbid[ct->nforbid][1]=m_base2;
			ct->AddForbiddenPair(i,m_base2);
			m_base2--;
			//ct->nforbid++;


		}
		else {
			//ct->npair++;
			//ct->pair[ct->npair][0]=i;
			//ct->pair[ct->npair][1]=m_base2;
			ct->AddPair(i,m_base2);
			m_base2--;
		}

	}
}

BOOL CPairDialog::OnInitDialog( ) {

	if (isforbid) {

		SetWindowText("Prohibit Base Pairs");
	}

	return CDialog::OnInitDialog();
}


