#if !defined(AFX_SEQUENCEEDITOR_H__B61902E1_BB2A_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_SEQUENCEEDITOR_H__B61902E1_BB2A_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// SequenceEditor.h : header file
//


#include "sequenceedit.h"
#include "Sequencedoc.h"

/////////////////////////////////////////////////////////////////////////////
// CSequenceEditor form view

#ifndef __AFXEXT_H__
#include <afxext.h>
#endif

struct CReadObject {

	CFormView *parent;
	char *texttoread;
	bool cancel;
	int start,stop,current;
	char *datapath;

};


class CSequenceEditor : public CFormView
{
protected:
	CSequenceEditor();           // protected constructor used by dynamic creation

	DECLARE_DYNCREATE(CSequenceEditor)

// Form Data
public:
	//{{AFX_DATA(CSequenceEditor)
	enum { IDD = IDD_SEQUENCEEDITOR_FORM };
	CString	m_Title;
	//}}AFX_DATA
	void OnInitialUpdate();
	CSequenceEdit *SequenceWindow;
	CEdit *CommentWindow;
	CFont *KpFont;
	CReadObject readobject;
	
	void Fold(bool isRNA);
	LRESULT OnReadDone(WPARAM wParam, LPARAM lParam); 
	bool ShouldSave();
	bool reading;
	LRESULT OnHighLight(WPARAM i, LPARAM j);
	bool readwhiletyping;


	CSequenceDoc* GetSeqDocument();
	
// Attributes
public:

// Operations
public:
	void Format();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CSequenceEditor)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	virtual ~CSequenceEditor();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
	//{{AFX_MSG(CSequenceEditor)
	afx_msg void OnFileSavesequence();
	afx_msg void OnFileSavesequenceas();
	afx_msg void OnChange();
	afx_msg void OnDestroy();
	afx_msg void OnFoldasdna();
	afx_msg void OnFoldasrna();
	afx_msg void OnReadback();
	afx_msg void OnReadReadwhiletyping();
	afx_msg void OnCopy();
	afx_msg void OnPaste();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_SEQUENCEEDITOR_H__B61902E1_BB2A_11D4_9F32_00C0F02A5F5D__INCLUDED_)
