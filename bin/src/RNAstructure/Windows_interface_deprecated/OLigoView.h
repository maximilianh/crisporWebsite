#if !defined(AFX_OLIGOVIEW_H__6030C308_0716_46D1_8C7D_A6400FFD510F__INCLUDED_)
#define AFX_OLIGOVIEW_H__6030C308_0716_46D1_8C7D_A6400FFD510F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// OLigoView.h : header file
//
#include "oligodoc.h"


/////////////////////////////////////////////////////////////////////////////
// COLigoView view

class COligoView : public CView
{
protected:
	COligoView();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(COligoView)

// Attributes
public:

	int current,column2,seqtextwidth,sequenceline,oligoline;
	bool overallandduplex,overall,duplex,target,intermolecular,intramolecular;
	CFont *symbol,*font,*subscript,*sequence;
	CButton *go,*left,*right,*tenleft,*tenright;
	int graphmin,graphmax;

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(COLigoView)
	protected:
	virtual void OnDraw(CDC* pDC);      // overridden to draw this view
	virtual void OnInitialUpdate();     // first time after construct
	//}}AFX_VIRTUAL

	
	COligoDoc *GetDocument();
	int TextOut(CString text,CDC* pDC, int x, int y, bool draw);
	int DGOut(CDC* pDC, int x, int y, bool draw) ;
	int fontsize,sequencesize;

	void OnLeft();
	void OnRight();
	void OnGo();
	void OnTenLeft();
	void OnTenRight();
	void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags); 
	void placecheck();
	int findheight(int graphzero,int bar,int graphtop, int graphbottom);
	CBrush poverallBrush,*potherBrush,predbrush;

	afx_msg void OnLButtonDown( UINT nFlags, CPoint point );

	void DrawBar(CDC *pDC, int i, int array, int x, int width, int graphzero, int graphtop, int graphbottom, CBrush *pBrush);



// Implementation
protected:
	virtual ~COligoView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
	//{{AFX_MSG(COLigoView)
		afx_msg void OnGraphFreeenergyOverall();
		afx_msg void OnGraphFreeenergyOverallandduplex();
	afx_msg void OnGraphFreeenergyBrokentargetstructure();
	afx_msg void OnGraphFreeenergyDuplex();
	afx_msg void OnGraphFreeenergyOligomerbimolecular();
	afx_msg void OnGraphFreeenergyOligomerunimolecular();
	afx_msg void OnLButtonDblClk(UINT nFlags, CPoint point);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_OLIGOVIEW_H__6030C308_0716_46D1_8C7D_A6400FFD510F__INCLUDED_)
