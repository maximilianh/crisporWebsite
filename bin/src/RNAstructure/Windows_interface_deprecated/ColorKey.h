#pragma once

#include "DrawView.h"
// CColorKey view

class CColorKey : public CView
{
	DECLARE_DYNCREATE(CColorKey)

protected:
	CColorKey();           // protected constructor used by dynamic creation
	virtual ~CColorKey();
	CFont *Font;

public:
	virtual void OnDraw(CDC* pDC);      // overridden to draw this view
	COLORREF *crColor; 
	void OnInitialUpdate();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:
	DECLARE_MESSAGE_MAP()
};


