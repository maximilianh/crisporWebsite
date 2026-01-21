// ColorKeyDoc.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "ColorKeyDoc.h"


// CColorKeyDoc

IMPLEMENT_DYNCREATE(CColorKeyDoc, CDocument)

CColorKeyDoc::CColorKeyDoc()
{
}

BOOL CColorKeyDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;
	return TRUE;
}

CColorKeyDoc::~CColorKeyDoc()
{
}


BEGIN_MESSAGE_MAP(CColorKeyDoc, CDocument)
END_MESSAGE_MAP()



// CColorKeyDoc diagnostics

#ifdef _DEBUG
void CColorKeyDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CColorKeyDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG


// CColorKeyDoc serialization

void CColorKeyDoc::Serialize(CArchive& ar)
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


// CColorKeyDoc commands
