// ColorKey.cpp : implementation file
//  Draw a key for color annotation of drawing

#include "stdafx.h"
#include "RNAstructure.h"
#include "ColorKey.h"

#define KEYFONTHEIGHT 18 //font size for the key

// CColorKey

IMPLEMENT_DYNCREATE(CColorKey, CView)

CColorKey::CColorKey()
{
	crColor = new COLORREF[9];
	crColor[0]=RGB(0,0,0);
		crColor[1]=RGB(colorintensity,0,0);
		crColor[2]=RGB(colorintensity,colorintensity/2,0);
		crColor[3]=RGB(4*colorintensity/5,3*colorintensity/4,0);
		crColor[4]=RGB(colorintensity/2,colorintensity/2,0);
		crColor[5]=RGB(0,colorintensity,0);
		crColor[6]=RGB(0,colorintensity,colorintensity);
		crColor[7]=RGB(0,0,colorintensity);
		crColor[8]=RGB(colorintensity,0,colorintensity);


	LOGFONT pLogFont;
	
	pLogFont.lfHeight = KEYFONTHEIGHT;
	pLogFont.lfWidth = 0;
	pLogFont.lfEscapement = 0;
	pLogFont.lfOrientation = 0;
	pLogFont.lfWeight = FW_BOLD;
	pLogFont.lfItalic = 0;
	pLogFont.lfUnderline = 0;
	pLogFont.lfStrikeOut = 0;
	pLogFont.lfCharSet = ANSI_CHARSET;
	pLogFont.lfOutPrecision = OUT_DEFAULT_PRECIS;
	pLogFont.lfClipPrecision = CLIP_DEFAULT_PRECIS;
	pLogFont.lfQuality = PROOF_QUALITY;
	pLogFont.lfPitchAndFamily = VARIABLE_PITCH|FF_MODERN;
	strcpy(pLogFont.lfFaceName,"Courier New");

	Font= new CFont;
	Font->CreateFontIndirect(&pLogFont);
}

CColorKey::~CColorKey()
{

	delete[] crColor;
	delete Font;
}

BEGIN_MESSAGE_MAP(CColorKey, CView)
END_MESSAGE_MAP()

void CColorKey::OnInitialUpdate() {
	//ResizeParentToFit();
	CView::OnInitialUpdate();
}
// CColorKey drawing

void CColorKey::OnDraw(CDC* pDC)
{
	int y;
	CDocument* pDoc = GetDocument();
	

	//Just write a key for color annotation:
	y = KEYFONTHEIGHT/2;

	//select the font
	CFont *oldFont = pDC->SelectObject(Font);

	COLORREF startcolor = SetTextColor(pDC->m_hDC ,crColor[0]);
	
	SetTextColor(pDC->m_hDC ,crColor[1]);
	pDC->TextOut(0, y, "Probability >= 99%");
	y+=3*KEYFONTHEIGHT/2;

	SetTextColor(pDC->m_hDC ,crColor[2]);
	pDC->TextOut(0, y, "99% > Probability >= 95%");
	y+=3*KEYFONTHEIGHT/2;

	SetTextColor(pDC->m_hDC ,crColor[3]);
	pDC->TextOut(0, y, "95% > Probability >= 90%");
	y+=3*KEYFONTHEIGHT/2;

	SetTextColor(pDC->m_hDC ,crColor[4]);
	pDC->TextOut(0, y, "90% > Probability >= 80%");
	y+=3*KEYFONTHEIGHT/2;

	SetTextColor(pDC->m_hDC ,crColor[5]);
	pDC->TextOut(0, y, "80% > Probability >= 70%");
	y+=3*KEYFONTHEIGHT/2;

	SetTextColor(pDC->m_hDC ,crColor[6]);
	pDC->TextOut(0, y, "70% > Probability >= 60%");
	y+=3*KEYFONTHEIGHT/2;

	SetTextColor(pDC->m_hDC ,crColor[7]);
	pDC->TextOut(0, y, "60% > Probability > 50%");
	y+=3*KEYFONTHEIGHT/2;

	SetTextColor(pDC->m_hDC ,crColor[8]);
	pDC->TextOut(0, y, "50% >= Probability");
	y+=3*KEYFONTHEIGHT/2;

	
	


	//Restore the old font
	pDC->SelectObject(oldFont);
	SetTextColor(pDC->m_hDC ,startcolor);

	

}


// CColorKey diagnostics

#ifdef _DEBUG
void CColorKey::AssertValid() const
{
	CView::AssertValid();
}

void CColorKey::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}
#endif //_DEBUG


// CColorKey message handlers
