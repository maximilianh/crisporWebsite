// RNAstructureDoc.cpp : implementation of the CRNAstructureDoc class
//

#include "stdafx.h"
#include "RNAstructure.h"

#include "RNAstructureDoc.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CRNAstructureDoc

IMPLEMENT_DYNCREATE(CRNAstructureDoc, CDocument)

BEGIN_MESSAGE_MAP(CRNAstructureDoc, CDocument)
	//{{AFX_MSG_MAP(CRNAstructureDoc)
		// NOTE - the ClassWizard will add and remove mapping macros here.
		//    DO NOT EDIT what you see in these blocks of generated code!
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CRNAstructureDoc construction/destruction

CRNAstructureDoc::CRNAstructureDoc()
{
	// TODO: add one-time construction code here

}

CRNAstructureDoc::~CRNAstructureDoc()
{
}

BOOL CRNAstructureDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;

	// TODO: add reinitialization code here
	// (SDI documents will reuse this document)

	return TRUE;
}



/////////////////////////////////////////////////////////////////////////////
// CRNAstructureDoc serialization

void CRNAstructureDoc::Serialize(CArchive& ar)
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
// CRNAstructureDoc diagnostics

#ifdef _DEBUG
void CRNAstructureDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CRNAstructureDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CRNAstructureDoc commands
