// OligoDoc.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "OligoDoc.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// COligoDoc

IMPLEMENT_DYNCREATE(COligoDoc, CDocument)

COligoDoc::COligoDoc()
{
}


COligoDoc::COligoDoc(COligoObject *Oligoobject,CMainFrame* pMainFrame, CWinApp *PParent) {


	oligoobject = Oligoobject;
	pParent = PParent;
	pmainframe = pMainFrame;
	
	


}


BOOL COligoDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;
	return TRUE;
}

COligoDoc::~COligoDoc()
{
}


BEGIN_MESSAGE_MAP(COligoDoc, CDocument)
	//{{AFX_MSG_MAP(COligoDoc)
		// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// COligoDoc diagnostics

#ifdef _DEBUG
void COligoDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void COligoDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// COligoDoc serialization

void COligoDoc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}
}

/////////////////////////////////////////////////////////////////////////////
// COligoDoc commands
