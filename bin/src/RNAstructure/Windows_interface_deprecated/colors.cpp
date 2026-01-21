// Number.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "colors.h"



#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CNumber dialog


//max is the largest number structure that can be drawn
CColors::CColors(PlotDoc *doc, CWnd* pParent /*=NULL*/)
	: CDialog(CColors::IDD, pParent)
{
	//{{AFX_DATA_INIT(CColors)
	m_number = doc->colors;
	//}}AFX_DATA_INIT

	Parent = pParent;
	Doc=doc;
}


void CColors::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CColors)
	DDX_Text(pDX, IDC_NUMBER, m_number);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CColors, CDialog)
	//{{AFX_MSG_MAP(CColors)
		// NOTE: the ClassWizard will add message map macros here
		ON_NOTIFY(UDN_DELTAPOS, IDC_SPIN, OnDeltaposSpin)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CNumber message handlers
BOOL CColors::OnInitDialog() {
	BOOL Return = CDialog::OnInitDialog();

	//Create the Spin Control
	CWnd* pWnd = GetDlgItem( IDC_SPIN_POS );
	CRect rect;
	pWnd->GetWindowRect( &rect );
	ScreenToClient( &rect );

	m_Spin.Create( WS_VISIBLE|WS_CHILD/*|dwStyles*/, rect, this, IDC_SPIN );
	m_Spin.SetRange( 1, 17 );  // Sends UDM_SETRANGE
	m_Spin.SetPos(m_number);


	return Return;
}

void CColors::OnDeltaposSpin(NMHDR* pNMHDR, LRESULT* pResult) 
{
	NM_UPDOWN* pNMUpDown = (NM_UPDOWN*)pNMHDR;
	// TODO: Add your control notification handler code here
	
	//increase the increment to a factor of ten
	
	m_number = m_number + 2*pNMUpDown->iDelta;
	
	if (m_number>2&&m_number<maxcolors) { 
		m_Spin.SetPos(m_number);
		UpdateData(FALSE);
		Doc->colors=m_number;
		Doc->colorranges();
		Parent->Invalidate();
		

	}
	else m_number = m_number - 2*pNMUpDown->iDelta; 


	*pResult = 0;
}