#if !defined(AFX_BIFOLDDOC_H__8AFD9841_D44E_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_BIFOLDDOC_H__8AFD9841_D44E_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// BiFoldDoc.h : header file
//

#include "../src/rna_library.h"
/////////////////////////////////////////////////////////////////////////////
// CBiFoldDoc document

class CBiFoldDoc : public CCTDoc
{
protected:
	//CBiFoldDoc();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CBiFoldDoc)

// Attributes
public:
	CBiFoldDoc(bool isRNA=true, char *datapath=NULL, CWinApp *pMainFrame1=NULL, char *startpath1=NULL, bool fold = true, bool *savefile = NULL);
	datatable data;
	bool OK;
	CWinApp *pMainFrame;
	char *startpath,*Datapath;//keep track of the path that the user last used to open a file
							//this is saved and restored from the registry
	CFrameWnd* Frame;
	CMainFrame *menuframe;
	bool ISRNA,*savefile;
	float T;//the temperature in K
	int newtemp();//change the data files for the new temperature
// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CBiFoldDoc)
	public:
	virtual void Serialize(CArchive& ar);   // overridden for document i/o
	protected:
	virtual BOOL OnNewDocument();
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CBiFoldDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
protected:
	//{{AFX_MSG(CBiFoldDoc)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_BIFOLDDOC_H__8AFD9841_D44E_11D4_9F32_00C0F02A5F5D__INCLUDED_)
