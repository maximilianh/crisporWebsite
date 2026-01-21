// GoDialog.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "GoDialog.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CGoDialog dialog


CGoDialog::CGoDialog(int *Current, COligoObject *Oligoobject, CWnd* pParent /*=NULL*/)
	: CDialog(CGoDialog::IDD, pParent)
{
	//{{AFX_DATA_INIT(CGoDialog)
	m_current = *Current;
	//}}AFX_DATA_INIT
	current = Current;
	oligoobject = Oligoobject;
	parent = pParent;
}


void CGoDialog::DoDataExchange(CDataExchange* pDX)
{
	
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CGoDialog)
	DDX_Text(pDX, IDC_NUMBER, m_current);
	DDV_MinMaxInt(pDX, m_current, 1, 100000);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CGoDialog, CDialog)
	ON_NOTIFY(UDN_DELTAPOS, IDC_SPIN, OnDeltaposSpin)
	//{{AFX_MSG_MAP(CGoDialog)
	ON_BN_CLICKED(IDMOSTSTABLE, OnMoststable)
	ON_EN_CHANGE(IDC_NUMBER, OnChangeNumber)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CGoDialog message handlers

void CGoDialog::OnOK() 
{
	UpdateData(TRUE);
	*current = m_current;

	parent->InvalidateRect(NULL);
	parent->UpdateWindow();
	
	CDialog::OnOK();
}

void CGoDialog::OnMoststable() 
{
	//find the most stable
	int i,ms;

	ms =1;
	for (i=2;i<=oligoobject->ct.GetSequenceLength()-oligoobject->length+1;i++) {
		if (oligoobject->table[i][0]<oligoobject->table[ms][0]) ms = i;

	}
	m_current = ms;
	UpdateData(FALSE);

	
}

void CGoDialog::OnChangeNumber() 
{
	// TODO: If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.
	
	// TODO: Add your control notification handler code here
	
}

BOOL CGoDialog::OnInitDialog( ) {

	BOOL Return = CDialog::OnInitDialog();

	//Create the Spin Control
	CWnd* pWnd = GetDlgItem( IDC_SPIN_POS );
	CRect rect;
	pWnd->GetWindowRect( &rect );
	ScreenToClient( &rect );

	m_Spin.Create( WS_VISIBLE|WS_CHILD/*|dwStyles*/, rect, this, IDC_SPIN );
	m_Spin.SetRange( 1, 10000 );  // Sends UDM_SETRANGE
	m_Spin.SetPos(m_current);


	return Return;
}

void CGoDialog::OnDeltaposSpin(NMHDR* pNMHDR, LRESULT* pResult) 
{
	NM_UPDOWN* pNMUpDown = (NM_UPDOWN*)pNMHDR;
	// TODO: Add your control notification handler code here
	
	//increase the increment to a factor of ten
	
	m_current = m_current + pNMUpDown->iDelta;
	
	if (m_current<=oligoobject->ct.GetSequenceLength()-oligoobject->length+1&&m_current>0) { 
		m_Spin.SetPos(m_current);
		UpdateData(FALSE);
		

	}
	else m_current = m_current - pNMUpDown->iDelta; 


	*pResult = 0;
}