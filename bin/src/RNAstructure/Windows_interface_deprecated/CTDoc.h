#if !defined(AFX_CTDOC_H__56CD8612_8D8A_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_CTDOC_H__56CD8612_8D8A_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// CTDoc.h : header file
//
#include "../src/structure.h"
/////////////////////////////////////////////////////////////////////////////
// CCTDoc document

class CCTDoc : public CDocument
{
protected:
	//CCTDoc();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CCTDoc)

// Attributes
public:
	CCTDoc();
	structure ct;
	int OpenSequence(char *filename);
	int integer(char nuc);


// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CCTDoc)
	public:
	virtual void Serialize(CArchive& ar);   // overridden for document i/o
	protected:
	virtual BOOL OnNewDocument();
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CCTDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
protected:
	//{{AFX_MSG(CCTDoc)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_CTDOC_H__56CD8612_8D8A_11D4_9F32_00C0F02A5F5D__INCLUDED_)
