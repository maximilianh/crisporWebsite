#if !defined(AFX_ZOOM_H__42541A61_B04F_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_ZOOM_H__42541A61_B04F_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// Zoom.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CZoom dialog

class CZoom : public CDialog
{
// Construction
public:
	CZoom(int *zoom, CWnd* pParent = NULL);   // standard constructor
	void OnOK();
	int *Zoom;
	CWnd* Parent;
	BOOL OnInitDialog();
	CSpinButtonCtrl m_Spin;

	

// Dialog Data
	//{{AFX_DATA(CZoom)
	enum { IDD = IDD_ZOOM };
	int		m_zoom;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CZoom)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CZoom)
	afx_msg void OnDeltaposSpin(NMHDR* pNMHDR, LRESULT* pResult);
	
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_ZOOM_H__42541A61_B04F_11D4_9F32_00C0F02A5F5D__INCLUDED_)
