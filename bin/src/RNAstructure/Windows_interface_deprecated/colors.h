#if !defined(AFX_NUMBER_H__DE6B7320_B3C5_11D4_9F32_00C0F02A5F5D__INCLUDED_2)
#define AFX_NUMBER_H__DE6B7320_B3C5_11D4_9F32_00C0F02A5F5D__INCLUDED_2

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// Number.h : header file
//

#include "plotdoc.h"


/////////////////////////////////////////////////////////////////////////////
// CNumber dialog

class CColors : public CDialog
{
// Construction
public:
	CColors(PlotDoc *doc, CWnd* pParent = NULL);   // standard constructor
	
	BOOL OnInitDialog();
	CSpinButtonCtrl m_Spin;
	PlotDoc *Doc;
	CWnd *Parent;
	
// Dialog Data
	//{{AFX_DATA(CColors)
	enum { IDD = IDD_COLORS };
	int		m_number;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CColors)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CColors)
		// NOTE: the ClassWizard will add member functions here
		afx_msg void OnDeltaposSpin(NMHDR* pNMHDR, LRESULT* pResult);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_NUMBER_H__DE6B7320_B3C5_11D4_9F32_00C0F02A5F5D__INCLUDED_2)
