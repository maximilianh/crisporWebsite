#if !defined(AFX_NUMBER_H__DE6B7320_B3C5_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_NUMBER_H__DE6B7320_B3C5_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// Number.h : header file
//
#include "drawdoc.h"


/////////////////////////////////////////////////////////////////////////////
// CNumber dialog

class CNumber : public CDialog
{
// Construction
public:
	CNumber(CDrawDoc *doc, CWnd* pParent = NULL);   // standard constructor
	
	BOOL OnInitDialog();
	CSpinButtonCtrl m_Spin;
	CDrawDoc *Doc;
	CWnd *Parent;
	
// Dialog Data
	//{{AFX_DATA(CNumber)
	enum { IDD = IDD_STRUCTURENUMBER };
	int		m_number;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CNumber)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	virtual void PostNcDestroy();
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CNumber)
		// NOTE: the ClassWizard will add member functions here
		afx_msg void OnDeltaposSpin(NMHDR* pNMHDR, LRESULT* pResult);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnEnChangeNumber();
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_NUMBER_H__DE6B7320_B3C5_11D4_9F32_00C0F02A5F5D__INCLUDED_)
