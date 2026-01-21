#if !defined(AFX_SEQUENCEDOC2_H__C0975BA1_BF2F_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_SEQUENCEDOC2_H__C0975BA1_BF2F_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// SequenceDoc2.h : header file
//
#include <fstream>
using namespace std;

/////////////////////////////////////////////////////////////////////////////
// CSequenceDoc2 document

class CSequenceDoc : public CDocument
{
protected:
	CSequenceDoc();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CSequenceDoc)

// Attributes
public:
	CSequenceDoc(CMainFrame* Frame, char *Startpath, char *Datapath,CWinApp *MainFrame, char *Filename=NULL);
	char *startpath, *filename,*datapath;
	bool filenamedefined;
	void allocatefilename(char *name, int length);
	CString title,comment,sequence;
	int OpenSequence(char *filename);
	int check;
	CWinApp *pMainFrame;
	bool changed;
	CMainFrame *pFrame;
	
	
// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CSequenceDoc2)
	public:
	virtual void Serialize(CArchive& ar);   // overridden for document i/o
	protected:
	virtual BOOL OnNewDocument();
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CSequenceDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
protected:
	//{{AFX_MSG(CSequenceDoc2)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_SEQUENCEDOC2_H__C0975BA1_BF2F_11D4_9F32_00C0F02A5F5D__INCLUDED_)
