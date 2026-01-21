#if !defined(AFX_SEQUENCEEDIT_H__8924A181_BB3F_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_SEQUENCEEDIT_H__8924A181_BB3F_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// SequenceEdit.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CSequenceEdit view

class CSequenceEdit : public CEdit
{
protected:
	//CSequenceEdit();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CSequenceEdit)

// Attributes
public:
	CSequenceEdit(bool *Readwhiletyping,char *datapath);
	CSequenceEdit();
	bool *read;
	char a[maxfil],c[maxfil],g[maxfil],u[maxfil],x[maxfil],t[maxfil];
	int foo;
	void OnSysChar(UINT nChar, UINT nRepCnt, UINT nFlags);
	void format();
// Operations
public:
	~CSequenceEdit();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CSequenceEdit)
	protected:
	
	//}}AFX_VIRTUAL

// Implementation
protected:
	//virtual ~CSequenceEdit();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
protected:
	//{{AFX_MSG(CSequenceEdit)
	afx_msg void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnChar(UINT nChar, UINT nRepCnt, UINT nFlags);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_SEQUENCEEDIT_H__8924A181_BB3F_11D4_9F32_00C0F02A5F5D__INCLUDED_)
