// OligoScreenDoc.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "OligoScreenDoc.h"


// COligoScreenDoc

IMPLEMENT_DYNCREATE(COligoScreenDoc, CDocument)

COligoScreenDoc::COligoScreenDoc()
{
}

COligoScreenDoc::COligoScreenDoc(char *DataPath, char* StartPath)
{
	datapath = DataPath;
	startpath = StartPath;
	SetTitle("OligoScreen");
	
}

BOOL COligoScreenDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;
	return TRUE;
}

COligoScreenDoc::~COligoScreenDoc()
{
}


BEGIN_MESSAGE_MAP(COligoScreenDoc, CDocument)
END_MESSAGE_MAP()


// COligoScreenDoc diagnostics

#ifdef _DEBUG
void COligoScreenDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void COligoScreenDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG


// COligoScreenDoc serialization

void COligoScreenDoc::Serialize(CArchive& ar)
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


// COligoScreenDoc commands
