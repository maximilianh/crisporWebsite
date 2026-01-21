// MultilignDoc.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "TurboFoldDoc.h"


// CMultilignDoc

IMPLEMENT_DYNCREATE(CTurboFoldDoc, CDocument)

CTurboFoldDoc::CTurboFoldDoc(bool isRNA, char *datapath1, CWinApp *pMainFrame1,char *startpath1)
{
	EnableAutomation();
	startpath = startpath1;

	OK = true;
	pMainFrame = pMainFrame1;
	ISRNA = isRNA;
	datapath = datapath1;

	if (ISRNA) SetTitle("RNA TurboFold");
	else SetTitle("DNA TurboFold"); 
}

BOOL CTurboFoldDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;
	return TRUE;
}



CTurboFoldDoc::~CTurboFoldDoc()
{
}

void CTurboFoldDoc::OnFinalRelease()
{
	// When the last reference for an automation object is released
	// OnFinalRelease is called.  The base class will automatically
	// deletes the object.  Add additional cleanup required for your
	// object before calling the base class.

	CDocument::OnFinalRelease();
}


BEGIN_MESSAGE_MAP(CTurboFoldDoc, CDocument)
END_MESSAGE_MAP()

BEGIN_DISPATCH_MAP(CTurboFoldDoc, CDocument)
END_DISPATCH_MAP()

// Note: we add support for IID_IMultilignDoc to support typesafe binding
//  from VBA.  This IID must match the GUID that is attached to the 
//  dispinterface in the .IDL file.

// {1114663F-4F72-44C0-A535-3A89B251D758}
static const IID IID_ITurboFoldDoc =
{ 0x1114663F, 0x4F72, 0x44C0, { 0xA5, 0x35, 0x3A, 0x89, 0xB2, 0x51, 0xD7, 0x58 } };

BEGIN_INTERFACE_MAP(CTurboFoldDoc, CDocument)
	INTERFACE_PART(CTurboFoldDoc, IID_ITurboFoldDoc, Dispatch)
END_INTERFACE_MAP()


// CMultilignDoc diagnostics

#ifdef _DEBUG
void CTurboFoldDoc::AssertValid() const
{
	CDocument::AssertValid();
}

#ifndef _WIN32_WCE
void CTurboFoldDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif
#endif //_DEBUG

#ifndef _WIN32_WCE
// CMultilignDoc serialization

void CTurboFoldDoc::Serialize(CArchive& ar)
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
#endif


// CMultilignDoc commands
