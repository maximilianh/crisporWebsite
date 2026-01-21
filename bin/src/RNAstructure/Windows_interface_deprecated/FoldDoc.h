#if !defined(AFX_FoldDoc_H__56CD8612_8D8A_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_FoldDoc_H__56CD8612_8D8A_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// FoldDoc.h : header file
//
#include "CTDoc.h"

#include "../src/rna_library.h"
#include "MainFrm.h"


/////////////////////////////////////////////////////////////////////////////
// CFoldDoc document

class CFoldDoc : public CCTDoc
{
protected:
	//CFoldDoc();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CFoldDoc)

// Attributes
public:
	CFoldDoc(bool isRNA=true, char *datapath=NULL, CWinApp *pMainFrame1=NULL, char *startpath1=NULL, bool fold = true, bool *Savefile= NULL, bool PF=false);
	datatable data;
	
	bool OK;
	CWinApp *pMainFrame;
	char* startpath;//keep track of the path that the user last used to open a file
							//this is saved and restored from the registry
	CFrameWnd* Frame;
	
	CMainFrame *menuframe;
	bool ISRNA;
	int maximuminternalloopsize;//the maximum number of unpaired nucleotides in an internal loop
	
	
	bool filenamedefined;
	char *filename,*Datapath;
	void namefile(char *Filename);
	bool *savefile;
	float T;//the temperature in K
	int newtemp();//change the data files for the new temperature
// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CFoldDoc)
	public:
	virtual void Serialize(CArchive& ar);   // overridden for document i/o
	protected:
	virtual BOOL OnNewDocument();
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CFoldDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
protected:
	//{{AFX_MSG(CFoldDoc)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_FoldDoc_H__56CD8612_8D8A_11D4_9F32_00C0F02A5F5D__INCLUDED_)
