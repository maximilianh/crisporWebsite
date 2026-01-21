#if !defined(AFX_OLIGODOC_H__EA0A96E4_049A_4032_A20F_6D3F49C26336__INCLUDED_)
#define AFX_OLIGODOC_H__EA0A96E4_049A_4032_A20F_6D3F49C26336__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// OligoDoc.h : header file
//


#include "../src/structure.h"
#include "oligoobject.h"

/////////////////////////////////////////////////////////////////////////////
// COligoDoc document

class COligoDoc : public CDocument
{
protected:
	COligoDoc();           // protected constructor used by dynamic creation
	
	DECLARE_DYNCREATE(COligoDoc)

// Attributes
public:

	COligoDoc(COligoObject *Oligoobject,CMainFrame* pMainFrame,CWinApp *PParent);

	
	CWinApp* pParent;

	CFrameWnd* Frame;
	COligoObject *oligoobject;
	CMainFrame* pmainframe;
	
	

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(COligoDoc)
	public:
	virtual void Serialize(CArchive& ar);   // overridden for document i/o
	protected:
	virtual BOOL OnNewDocument();
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~COligoDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
protected:
	//{{AFX_MSG(COligoDoc)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_OLIGODOC_H__EA0A96E4_049A_4032_A20F_6D3F49C26336__INCLUDED_)
