// ReFoldDoc.h: interface for the CReFoldDoc class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_REFOLDDOC_H__A66C28DB_44FC_4B11_AD40_B521DF66DEB6__INCLUDED_)
#define AFX_REFOLDDOC_H__A66C28DB_44FC_4B11_AD40_B521DF66DEB6__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "CTDoc.h"

class CReFoldDoc : public CCTDoc  
{
public:
	CReFoldDoc();
	virtual ~CReFoldDoc();
	char* startingpath; //startpath used in CRefoldView
	CFrameWnd* Frame;
	CWinApp *pMainFrame;

protected:
	
	DECLARE_DYNCREATE(CReFoldDoc)

};

#endif // !defined(AFX_REFOLDDOC_H__A66C28DB_44FC_4B11_AD40_B521DF66DEB6__INCLUDED_)
