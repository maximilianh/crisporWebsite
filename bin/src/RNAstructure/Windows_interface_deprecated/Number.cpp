// Number.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "Number.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CNumber dialog


//max is the largest number structure that can be drawn
CNumber::CNumber(CDrawDoc *doc, CWnd* pParent /*=NULL*/)
	: CDialog(CNumber::IDD, pParent)
{
	//{{AFX_DATA_INIT(CNumber)
	m_number = doc->StructureNumber;
	//}}AFX_DATA_INIT

	Parent = pParent;
	Doc=doc;
}


void CNumber::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CNumber)
	DDX_Text(pDX, IDC_NUMBER, m_number);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CNumber, CDialog)
	//{{AFX_MSG_MAP(CNumber)
		// NOTE: the ClassWizard will add message map macros here
		ON_NOTIFY(UDN_DELTAPOS, IDC_SPIN, OnDeltaposSpin)
	//}}AFX_MSG_MAP
	ON_EN_CHANGE(IDC_NUMBER, OnEnChangeNumber)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CNumber message handlers
BOOL CNumber::OnInitDialog() {
	BOOL Return = CDialog::OnInitDialog();

	//Create the Spin Control
	CWnd* pWnd = GetDlgItem( IDC_SPIN_POS );
	CRect rect;
	pWnd->GetWindowRect( &rect );
	ScreenToClient( &rect );

	m_Spin.Create( WS_VISIBLE|WS_CHILD/*|dwStyles*/, rect, this, IDC_SPIN );
	m_Spin.SetRange( 1, Doc->ct.GetNumberofStructures() );  // Sends UDM_SETRANGE
	m_Spin.SetPos(m_number);


	return Return;
}

void CNumber::OnDeltaposSpin(NMHDR* pNMHDR, LRESULT* pResult) 
{
	NM_UPDOWN* pNMUpDown = (NM_UPDOWN*)pNMHDR;
	// TODO: Add your control notification handler code here
	
	//increase the increment to a factor of ten
	
	m_number = m_number + pNMUpDown->iDelta;
	
	if (m_number<=Doc->ct.GetNumberofStructures()&&m_number>0) { 
		m_Spin.SetPos(m_number);
		UpdateData(FALSE);
		Doc->NewStructure(m_number);
		Parent->Invalidate();
		

	}
	else m_number = m_number - pNMUpDown->iDelta; 


	*pResult = 0;
}
void CNumber::OnEnChangeNumber()
{
	//the number has been changed
	UpdateData(TRUE);
	if (m_number<=Doc->ct.GetNumberofStructures()&&m_number>0) { 
		m_Spin.SetPos(m_number);
		
		Doc->NewStructure(m_number);
		Parent->Invalidate();

	}

}

void  CNumber::PostNcDestroy()
{
	delete this;
}

