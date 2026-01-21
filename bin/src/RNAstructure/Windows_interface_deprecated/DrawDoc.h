#if !defined(AFX_DRAWDOC_H__76570683_9982_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_DRAWDOC_H__76570683_9982_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// DrawDoc.h : header file
//
#include "ctdoc.h"
#include "../src/draw.h"
#include "../src/findpseudo.h"
#include "mainfrm.h"
#include "RNAstructure.h"


/////////////////////////////////////////////////////////////////////////////
// CDrawDoc document

class CDrawDoc : public CCTDoc
{
protected:
	
	DECLARE_DYNCREATE(CDrawDoc)
	CDrawDoc();

// Attributes
public:
	CDrawDoc(const char *ctfilename, bool *Clockwise, CMainFrame* pMainFrame, CRNAstructureApp *App);
	int FontSize;
	short StructureNumber;
	void Calculate();
	void NewStructure(int number);
	coordinates *out;
	int pseudo[maxpseudo];
	int npseudo;
	int height,width,xmax,ymax;
	bool nopair,*clockwise;
	CMainFrame* pmainframe;
	bool iscolorannotated,isSHAPEannotated,PK;
	short *color;
	void colorannotate();
	void OutPutHelices(char *filename);
	CRNAstructureApp *app;

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDrawDoc)
	public:
	virtual void Serialize(CArchive& ar);   // overridden for document i/o
	protected:
	virtual BOOL OnNewDocument();
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CDrawDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
protected:
	//{{AFX_MSG(CDrawDoc)
		// NOTE - the ClassWizard will add and remove member functions here.
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_DRAWDOC_H__76570683_9982_11D4_9F32_00C0F02A5F5D__INCLUDED_)
