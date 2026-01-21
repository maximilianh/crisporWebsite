#pragma once

#include "DrawView.h"
// CColorKey view

class CSHAPEKey : public CView
{
	DECLARE_DYNCREATE(CSHAPEKey)

protected:
	CSHAPEKey();           // protected constructor used by dynamic creation
	virtual ~CSHAPEKey();
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


