#if !defined(AFX_DRAWVIEW_H__9C9028E1_9A0E_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_DRAWVIEW_H__9C9028E1_9A0E_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// DrawView2.h : header file
//

#include "DrawDoc.h"
#include "number.h"
#include "../src/pfunction.h"
#include "../src/phmm/utils/xmath/log/xlog_math.h"


#define colorintensity 255 //intenity of colors used in color annotation

/////////////////////////////////////////////////////////////////////////////
// CDrawView2 view

class CDrawView : public CScrollView
{
protected:
	CDrawView();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CDrawView)

// Attributes
public:
	void NewFont();
	void NewSize(int Size);
	void NewStructure(int number);
	CDrawDoc *pDoc;
	int zoom;
	void Cm_UP();
	void Cm_DOWN();
	void Cm_LARGER();
	void Cm_SMALLER();
	bool resize;
	BOOL OnPreparePrinting(CPrintInfo *pinfo);
	void OnBeginPrinting(CDC *pDC, CPrintInfo *pinfo);
	void OnEndPrinting(CDC *pDC, CPrintInfo *pInfo);
	void OnPrint(CDC *pDC, CPrintInfo *pInfo);
	void draw(CDC* pDC);
	COLORREF *crColor;
	pfunctionclass *v;//needed for probability color annotation
	PFPRECISION *w5;//needed for probability color annotation
	PFPRECISION scaling;//needed for probability color annotation
	pfdatatable *data;//needed for probability color annotation
	bool *lfce,*mod;//needed for probability color annotation
	forceclass *fce;//needed for probability color annotation

	void SetColors();
	
	
	CFont *Font;
// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDrawView)
	protected:
	virtual void OnDraw(CDC* pDC);      // overridden to draw this view
	virtual void OnInitialUpdate();     // first time after construct
	//}}AFX_VIRTUAL

// Implementation
protected:
	virtual ~CDrawView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
	//{{AFX_MSG(CDrawView)
	afx_msg void OnDrawZoom();
	afx_msg void OnDrawStructurenumber();
	afx_msg void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnDrawClockwise();
	afx_msg void OnFilePrint();
	afx_msg void OnFilePrintPreview();
	afx_msg void OnCpy();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnDrawAddcolorannotation();
	afx_msg void OnDrawShowcolorannotationkey();
	afx_msg void OnDrawExportstructuretotextfile();
	void DetermineColor();
	void DetermineSHAPEColor();
	afx_msg void OnDrawAddshapeannotation();
	afx_msg void OnDrawShowshapeannotationkey();
	afx_msg void OnDrawExporttodot();
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_DRAWVIEW_H__9C9028E1_9A0E_11D4_9F32_00C0F02A5F5D__INCLUDED_)
