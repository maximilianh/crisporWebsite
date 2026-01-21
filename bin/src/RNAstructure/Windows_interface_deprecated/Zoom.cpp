// Zoom.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "Zoom.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif



/////////////////////////////////////////////////////////////////////////////
// CZoom dialog


CZoom::CZoom(int *zoom, CWnd* pParent /*=NULL*/)
	: CDialog(CZoom::IDD, pParent)
{
	//{{AFX_DATA_INIT(CZoom)
	m_zoom = *zoom;
	//}}AFX_DATA_INIT
	Zoom = zoom;
	Parent = pParent;
}


void CZoom::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CZoom)
	DDX_Text(pDX, IDC_ZOOM, m_zoom);
	DDV_MinMaxInt(pDX, m_zoom, 1, 10000);
	//}}AFX_DATA_MAP
}

void CZoom::OnOK() 
{
	
	UpdateData(TRUE);
	*Zoom = m_zoom;
	Parent->Invalidate();
	CDialog::OnOK();
}


BEGIN_MESSAGE_MAP(CZoom, CDialog)
	//{{AFX_MSG_MAP(CZoom)
	ON_NOTIFY(UDN_DELTAPOS, IDC_SPIN, OnDeltaposSpin)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CZoom message handlers



BOOL CZoom::OnInitDialog() {
	BOOL Return = CDialog::OnInitDialog();

	//Create the Spin Control
	CWnd* pWnd = GetDlgItem( IDC_SPIN_POS );
	CRect rect;
	pWnd->GetWindowRect( &rect );
	ScreenToClient( &rect );

	m_Spin.Create( WS_VISIBLE|WS_CHILD/*|dwStyles*/, rect, this, IDC_SPIN );
	m_Spin.SetRange( 1, 10000 );  // Sends UDM_SETRANGE
	m_Spin.SetPos(m_zoom);


	return Return;
}

void CZoom::OnDeltaposSpin(NMHDR* pNMHDR, LRESULT* pResult) 
{
	NM_UPDOWN* pNMUpDown = (NM_UPDOWN*)pNMHDR;
	// TODO: Add your control notification handler code here
	
	//increase the increment to a factor of ten
	
	m_zoom = m_zoom + 5*pNMUpDown->iDelta;
	
	if (m_zoom<10000&&m_zoom>0) { 
		m_Spin.SetPos(m_zoom);
		UpdateData(FALSE);
		*Zoom = m_zoom;
		Parent->Invalidate();

	}
	else m_zoom = m_zoom - 5*pNMUpDown->iDelta; 


	*pResult = 0;
}