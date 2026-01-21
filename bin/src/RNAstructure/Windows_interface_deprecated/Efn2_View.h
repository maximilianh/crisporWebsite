#if !defined(AFX_EFN2_VIEW_H__E73BF322_CA05_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_EFN2_VIEW_H__E73BF322_CA05_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// Efn2_View.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CEfn2_View form view

#ifndef __AFXEXT_H__
#include <afxext.h>
#endif

class CEfn2_View : public CFormView
{
protected:
	CEfn2_View();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CEfn2_View)

// Form Data
public:
	//{{AFX_DATA(CEfn2_View)
	enum { IDD = IDD_EFN2_VIEW_FORM };
	CString	m_ctfilename;
	CString	m_outfilename;
	//}}AFX_DATA
	CFoldDoc *GetFoldDocument();
	void OnTemperature();

// Attributes
public:

// Operations
public:
	void OnInitialUpdate();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CEfn2_View)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	virtual ~CEfn2_View();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
	//{{AFX_MSG(CEfn2_View)
	afx_msg void OnCtfileEfn();
	afx_msg void OnOutfileEfn();
	afx_msg void OnStartefn();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:
	BOOL m_details;
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_EFN2_VIEW_H__E73BF322_CA05_11D4_9F32_00C0F02A5F5D__INCLUDED_)
