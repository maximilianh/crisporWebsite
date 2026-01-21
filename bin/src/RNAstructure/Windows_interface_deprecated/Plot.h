#pragma once
#include "plotdoc.h"


// CPlot view

class CPlot : public CScrollView
{
	DECLARE_DYNCREATE(CPlot)

protected:
	CPlot();           // protected constructor used by dynamic creation
	virtual ~CPlot();

public:
	virtual void OnDraw(CDC* pDC);      // overridden to draw this view
	PlotDoc *plotdoc;
	void OnInitialUpdate();
	
	int zoom;

	CSize *PageSize;
	CSize *WindowSize;
	bool firstdraw;
	int plotsize;
	CFont Font,Font2;
	afx_msg void OnLButtonDown( UINT nFlags, CPoint point ); 

	void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags);
	
	



#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnDrawZoom();
	afx_msg void OnDrawColors();
	afx_msg void OnDrawShowlegend();
	afx_msg void OnDrawPlotrange();
	void OnFilePrint(); 
	void OnFilePrintPreview();
	BOOL OnPreparePrinting(CPrintInfo* pInfo);
	void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	void OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/);
	void OnPrint( CDC *pDC, CPrintInfo *pInfo );
	void Draw(CDC* pDC);
	afx_msg void OnOutputplot();
};