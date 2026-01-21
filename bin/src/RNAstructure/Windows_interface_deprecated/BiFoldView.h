#if !defined(AFX_BIFOLDVIEW_H__8AFD9840_D44E_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_BIFOLDVIEW_H__8AFD9840_D44E_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// BiFoldView.h : header file
//
#include "bifolddoc.h"
#include "tprogressdialog.h"
#include "foldview.h"

/////////////////////////////////////////////////////////////////////////////
// CBiFoldView dialog

class CBiFoldView : public CFormView
{
// Construction
public:
	CBiFoldView(CWnd* pParent = NULL);   // standard constructor
	CBiFoldDoc *GetBiFoldDocument();
	void OnUpdate(CView*, LPARAM, CObject*);
	structure ct2,ct3;
	int dynwindow,dynnumber,dynpercent;
	TProgressDialog *progress;
	CFoldObject foldobject;
	LRESULT DoneFolding(WPARAM wParam, LPARAM lParam);
	//CFoldDoc *GetFoldDocument();
	void OnInitialUpdate(); 
	bool forbidunimolecular,started;
	void OnTemperature();
// Dialog Data
	//{{AFX_DATA(CBiFoldView)
	enum { IDD = IDD_BIFOLDVIEW_FORM };
	CString	m_ctname;
	CString	m_sequence1;
	CString	m_sequence2;
	short	m_window;
	short	m_percent;
	short	m_number;
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CBiFoldView)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	DECLARE_DYNCREATE(CBiFoldView)
	//CBiFoldView(); //protected constructor 

	// Generated message map functions
	//{{AFX_MSG(CBiFoldView)
	afx_msg void OnSequencebutton();
	afx_msg void OnCtbutton();
	afx_msg void OnSequence2button();
	afx_msg void OnStart();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:
	BOOL m_save;
	afx_msg void OnForceForbidunimolecularpairs();
	afx_msg void OnSetFocus(CWnd* pOldWnd);
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_BIFOLDVIEW_H__8AFD9840_D44E_11D4_9F32_00C0F02A5F5D__INCLUDED_)
