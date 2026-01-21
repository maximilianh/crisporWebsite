// ColorKey.cpp : implementation file
//  Draw a key for color annotation of drawing

#include "stdafx.h"
#include "RNAstructure.h"
#include "SHAPEColorKey.h"

#define KEYFONTHEIGHT 18 //font size for the key

// CColorKey

IMPLEMENT_DYNCREATE(CSHAPEKey, CView)

CSHAPEKey::CSHAPEKey()
{
	crColor = new COLORREF[9];

	//set the key to no data = grey; >.7 red; .3 to .7 orange; <.3 black

	crColor[0]=RGB(0,0,0);//black
	crColor[1]=RGB(colorintensity,0,0);//red
	crColor[2]=RGB(colorintensity,colorintensity/2,0);//orange
	crColor[3]=RGB(colorintensity*2/3,colorintensity*2/3,colorintensity*2/3);//grey	

		//retired colors
		//crColor[3]=RGB(4*colorintensity/5,3*colorintensity/4,0);
		//crColor[4]=RGB(colorintensity/2,colorintensity/2,0);
		//crColor[5]=RGB(0,colorintensity,0);
		//crColor[6]=RGB(0,colorintensity,colorintensity);
		//crColor[7]=RGB(0,0,colorintensity);
		//crColor[8]=RGB(colorintensity,0,colorintensity);


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

CSHAPEKey::~CSHAPEKey()
{

	delete[] crColor;
	delete Font;
}

BEGIN_MESSAGE_MAP(CSHAPEKey, CView)
END_MESSAGE_MAP()

void CSHAPEKey::OnInitialUpdate() {
	//ResizeParentToFit();
	CView::OnInitialUpdate();
}
// CColorKey drawing

void CSHAPEKey::OnDraw(CDC* pDC)
{
	int y;
	CDocument* pDoc = GetDocument();
	

	//Just write a key for color annotation:
	y = KEYFONTHEIGHT/2;

	//select the font
	CFont *oldFont = pDC->SelectObject(Font);

	COLORREF startcolor = SetTextColor(pDC->m_hDC ,crColor[3]);
	pDC->TextOut(0,y,"No Data");
	y+=3*KEYFONTHEIGHT/2;
	
	SetTextColor(pDC->m_hDC ,crColor[1]);
	pDC->TextOut(0, y, "SHAPE >= 0.85");
	y+=3*KEYFONTHEIGHT/2;

	SetTextColor(pDC->m_hDC ,crColor[2]);
	pDC->TextOut(0, y, "0.85 > SHAPE >= 0.4");
	y+=3*KEYFONTHEIGHT/2;

	//SetTextColor(pDC->m_hDC ,crColor[3]);
	//pDC->TextOut(0, y, "1.5 > SHAPE >= 1.0");
	//y+=3*KEYFONTHEIGHT/2;

	//SetTextColor(pDC->m_hDC ,crColor[4]);
	//pDC->TextOut(0, y, "90% > BP Probability >= 80%");
	//y+=3*KEYFONTHEIGHT/2;

	//SetTextColor(pDC->m_hDC ,crColor[5]);
	//pDC->TextOut(0, y, "1.0 > SHAPE >= 0.5");
	//y+=3*KEYFONTHEIGHT/2;

	//SetTextColor(pDC->m_hDC ,crColor[6]);
	//pDC->TextOut(0, y, "0.5 > SHAPE >= 0.0");
	//y+=3*KEYFONTHEIGHT/2;

	SetTextColor(pDC->m_hDC ,crColor[0]);
	pDC->TextOut(0, y, "0.4 > SHAPE ");
	y+=3*KEYFONTHEIGHT/2;

	//SetTextColor(pDC->m_hDC ,crColor[8]);
	//pDC->TextOut(0, y, "50% >= BP Probability");
	//y+=3*KEYFONTHEIGHT/2;

	
	


	//Restore the old font
	pDC->SelectObject(oldFont);
	SetTextColor(pDC->m_hDC ,startcolor);

	

}


// CColorKey diagnostics

#ifdef _DEBUG
void CSHAPEKey::AssertValid() const
{
	CView::AssertValid();
}

void CSHAPEKey::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}
#endif //_DEBUG


// CColorKey message handlers
