// Plot.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "Plot.h"
#include "plotrange.h"
#include "../src/algorithm.h"

#define topmargin 50
#define bottommargin 10
#define leftmargin 10
#define rightmargin 50
#define dotsize 10
#define dotarea 8
#define extlinewidth 6
#define intlinewidth 2
#define FontSize 32
#define scalespace 10
#define fontleftshift 6
#define lineext 8
#include "zoom.h"
#include "colors.h"
#include <math.h>



#define SCROLLCHANGE 5

using namespace std;

// CPlot

IMPLEMENT_DYNCREATE(CPlot, CScrollView)

CPlot::CPlot()
{

	firstdraw = true;
}

CPlot::~CPlot()
{
	
	delete PageSize;
	delete WindowSize;
}


void CPlot::OnInitialUpdate()
{
	
	CString title;
	
	
	CSize sizeTotal;
	

	CScrollView::OnInitialUpdate();

	
	
	CRect rect;
	CPoint topleft,bottomright;
	
	//CSize sizeTotal;
	plotdoc = ((PlotDoc*) GetDocument());


	
	//////////////////////////
	//make sure the menu item for show legend is correctly checked
	CMenu* menu = plotdoc->app->pMainFrame->GetMenu();
	if (plotdoc->showlegend) menu->CheckMenuItem(ID_DRAW_SHOWLEGEND,MF_CHECKED);
	else menu->CheckMenuItem(ID_DRAW_SHOWLEGEND,MF_UNCHECKED); 

		
	//zoom out so that the whole structure fits in the client area
	GetClientRect(rect);

	rect.bottom;
	if (rect.bottom<rect.right) {

		zoom = int (100.0*(((float) (rect.bottom)) / ((float) (plotdoc->ct.GetSequenceLength()*dotsize+topmargin+bottommargin))) );
		
	}
	else zoom = int (100.0*(((float) (rect.right)) / ((float) (plotdoc->ct.GetSequenceLength()*dotsize+topmargin+bottommargin))) );
	
	if (zoom>250) zoom =250;

	PageSize = new CSize(int ((float)(plotdoc->ct.GetSequenceLength()*dotsize+topmargin+bottommargin)*(((float) (zoom)/(float) (100)))),
   		int ((float)(plotdoc->ct.GetSequenceLength()*dotsize+topmargin+bottommargin)*(((float) (zoom)/(float) (100)))));
	WindowSize= new CSize(plotdoc->ct.GetSequenceLength()*dotsize+topmargin+bottommargin,plotdoc->ct.GetSequenceLength()*dotsize+topmargin+bottommargin);
	
	//set the scroll bar sizes:
	//sizeTotal.cx = pDoc->xmax*((long (zoom)/long (100)))+15;
	//sizeTotal.cy = pDoc->ymax*((long (zoom)/long (100)))+15;
	//SetScrollSizes(MM_TEXT, sizeTotal);
	sizeTotal.cx = PageSize->cx+15;
	sizeTotal.cy = PageSize->cy+15+(plotdoc->colors)*FontSize/2;
	SetScrollSizes(MM_TEXT, sizeTotal);

	plotsize = (plotdoc->ct.GetSequenceLength()+1)*dotsize;


	///////////
	//Select a font for the scale info
	LOGFONT pLogFont;
	
	pLogFont.lfHeight = FontSize;
	pLogFont.lfWidth = 0;
	pLogFont.lfEscapement = 0;
	pLogFont.lfOrientation = 0;
	pLogFont.lfWeight = FW_NORMAL;
	pLogFont.lfItalic = 0;
	pLogFont.lfUnderline = 0;
	pLogFont.lfStrikeOut = 0;
	pLogFont.lfCharSet = ANSI_CHARSET;
	pLogFont.lfOutPrecision = OUT_DEFAULT_PRECIS;
	pLogFont.lfClipPrecision = CLIP_DEFAULT_PRECIS;
	pLogFont.lfQuality = PROOF_QUALITY;
	pLogFont.lfPitchAndFamily = VARIABLE_PITCH|FF_MODERN;
	strcpy(pLogFont.lfFaceName,"Courier New");

	Font.CreateFontIndirect(&pLogFont);


	///////////////
	//Select a font for the legend
	
	pLogFont.lfHeight = FontSize/2;
	

	Font2.CreateFontIndirect(&pLogFont);
	
	
}

BEGIN_MESSAGE_MAP(CPlot, CScrollView)
	ON_COMMAND(ID_DRAW_ZOOM, OnDrawZoom)
	ON_WM_KEYDOWN()
	ON_WM_LBUTTONDOWN()
	ON_COMMAND(ID_DRAW_COLORS, OnDrawColors)
	ON_COMMAND(ID_DRAW_SHOWLEGEND, OnDrawShowlegend)
	ON_COMMAND(ID_DRAW_PLOTRANGE, OnDrawPlotrange)
	ON_COMMAND(ID_FILE_PRINT, OnFilePrint)
	ON_COMMAND(ID_FILE_PRINTPREVIEW, OnFilePrintPreview)
	ON_COMMAND(ID_OUTPUTPLOT, OnOutputplot)
END_MESSAGE_MAP()


// CPlot drawing

void CPlot::OnDraw(CDC* pDC)
{

	RECT rect;
	GetClientRect(&rect);
	
	
	int prevMode;
	CSize  oldVExt, oldWExt;
	CSize sizeTotal;


	//////////////////////////////////////////////
	////////Set the Viewport and Window sizes



	prevMode = pDC->SetMapMode(MM_ANISOTROPIC);

	PageSize = new CSize(int ((float)(plotdoc->ct.GetSequenceLength()*dotsize+topmargin+bottommargin)*(((float) (zoom)/(float) (100)))),
   		int ((float)(plotdoc->ct.GetSequenceLength()*dotsize+topmargin+bottommargin)*(((float) (zoom)/(float) (100)))));
	WindowSize= new CSize(plotdoc->ct.GetSequenceLength()*dotsize+topmargin+bottommargin,plotdoc->ct.GetSequenceLength()*dotsize+topmargin+bottommargin);

	sizeTotal.cx = PageSize->cx+15;
	sizeTotal.cy = PageSize->cy+15+(plotdoc->colors)*FontSize/2;
	SetScrollSizes(MM_TEXT, sizeTotal);

	oldVExt = pDC->SetViewportExt(PageSize->cx,PageSize->cy);

	oldWExt = pDC->SetWindowExt(WindowSize->cx,WindowSize->cy); 
	
	
	
	//shrink the view if the dot plot is smaller when first displayed 
	if (firstdraw) {
		ResizeParentToFit();
		firstdraw = false;
	}

	

	///////////////////////////////////////////////////////
	/////////Establish the total size of the plot
	/*int xplotsize,yplotsize;

	if (((rect.right-rightmargin-leftmargin)/(plotdoc->xstop-plotdoc->xstart))<((rect.bottom-topmargin-bottommargin)/(plotdoc->ystop-plotdoc->ystart))) {
		//the x dimension is limiting
		xplotsize = rect.right-rightmargin-leftmargin;
		yplotsize = (int)((float) xplotsize*(((float) (plotdoc->ystop-plotdoc->ystart))/((float) (plotdoc->xstop-plotdoc->xstart))));

	}
	else {
		//the y dimension is limiting
		yplotsize = rect.bottom-topmargin-bottommargin;
		xplotsize = (int) ((float) yplotsize*(((float) (plotdoc->xstop-plotdoc->xstart))/(float) (plotdoc->ystop-plotdoc->ystart)));

	}*/

	
	Draw(pDC);

	


	/////////////////////////////////////////////////
	//Restore Previous Window Settings
	pDC->SetViewportExt(oldVExt);
	pDC->SetWindowExt(oldWExt);
	pDC->SetMapMode(prevMode);
	
	
	
}

void CPlot::Draw(CDC* pDC) {
	COLORREF black = 0;
	CPen *pen;
	int i,j,k;

	/////////////////////////////////////////////////////////////////
	//Print a line around the plot (triangular) area:
		//use a extlinewidth pixel border:
	pen = new CPen;
	
	pen->CreatePen(PS_SOLID,extlinewidth,black);
	pDC->SelectObject(pen);

	pDC->MoveTo(leftmargin,topmargin);
	pDC->LineTo(plotsize+leftmargin,topmargin);
	pDC->LineTo(plotsize+leftmargin,plotsize+topmargin);
	pDC->LineTo(leftmargin,topmargin);
	delete pen;

	///////////////////////////////////////////////////////////////
	//Print lines every 10 nucs throughout plot
	pen = new CPen;
	pen->CreatePen(PS_SOLID,intlinewidth,black);
	pDC->SelectObject(pen);
	

	for (i=10;i<plotdoc->ct.GetSequenceLength();i+=10) {

		pDC->MoveTo(plotsize+leftmargin+extlinewidth+lineext/2,topmargin+dotsize*i+extlinewidth-dotsize/2);
		pDC->LineTo(leftmargin+i*dotsize-extlinewidth/2+dotsize/2,topmargin+dotsize*i+extlinewidth-dotsize/2);
		pDC->MoveTo(leftmargin+i*dotsize,topmargin+dotsize*i+extlinewidth-dotsize/2);
		pDC->LineTo(leftmargin+i*dotsize,topmargin-lineext);

	}

	delete pen;
	//////////////////////////////////////////////////////////////////
	//Add a scale to the plot:
	CFont *oldFont = pDC->SelectObject(&Font);
	char number[20];
	CString string;


	itoa(1,number,10);
	string = number;
	pDC->TextOut(leftmargin-fontleftshift, topmargin-FontSize-scalespace, string);
	pDC->TextOut(leftmargin+plotsize+extlinewidth+scalespace,topmargin+extlinewidth-FontSize/2,string);

	for (i=10;i<=plotdoc->ct.GetSequenceLength()-5;i+=10) {

		itoa(i,number,10);
		string=number;
		pDC->TextOut(leftmargin+(i-1)*dotsize-fontleftshift,topmargin-FontSize-scalespace, string);
		pDC->TextOut(leftmargin+plotsize+extlinewidth+scalespace,topmargin+extlinewidth+(i-1)*dotsize-FontSize/2,string);


	}

	itoa(plotdoc->ct.GetSequenceLength(),number,10);
	string = number;
	pDC->TextOut(leftmargin+(plotdoc->ct.GetSequenceLength()-1)*dotsize-fontleftshift, topmargin-FontSize-scalespace, string);
	pDC->TextOut(leftmargin+plotsize+extlinewidth+scalespace,topmargin+extlinewidth+(plotdoc->ct.GetSequenceLength()-1)*dotsize-FontSize/2,string);


	////////////////////////////////////////////////////////////////
	////Place Legend in lower left-hand corner:

	if (plotdoc->showlegend) {
		pDC->SelectObject(&Font2);

		for (i=0;i<plotdoc->colors;i++) {
			gcvt(plotdoc->colorrange[i],6,number);
			string = number;
			string+= " < ";
			string+=plotdoc->message;
			string+=" < ";
			gcvt(plotdoc->colorrange[i+1],6,number);
			string+=number;
			pDC->TextOut(leftmargin+dotsize,plotsize+topmargin+extlinewidth+(plotdoc->colors)*FontSize/2-(plotdoc->colors-i)*FontSize/2,string);
		
		}

		for (i=0;i<plotdoc->colors;i++) {
			pDC->FillSolidRect(leftmargin,
						plotsize+topmargin+extlinewidth+(plotdoc->colors)*FontSize/2-(plotdoc->colors-i)*FontSize/2,
						dotsize,dotsize,plotdoc->rgbcolors[i]);
		}
	}
	/////////////////////////////////////////////////////////////////
	/////add points:


	//Fill dots:
	
	for (i=plotdoc->xstart;i<=plotdoc->xstop;i++) {
		for (j=max(i+1,plotdoc->ystart);j<=plotdoc->ystop;j++) {
			for (k=1;k<=plotdoc->colors;k++) {
				if (plotdoc->arrayvalues[j][i]<=plotdoc->colorrange[k]&&plotdoc->arrayvalues[j][i]>=plotdoc->colorrange[k-1]) {
					pDC->FillSolidRect(leftmargin + extlinewidth + (j-plotdoc->xstart)*dotsize,
						topmargin+extlinewidth +(i-plotdoc->ystart)*dotsize,
						dotarea,dotarea,plotdoc->rgbcolors[k-1]);	
					break;
				}
			}
		}

	}

	pDC->SelectObject(oldFont);

}


// CPlot diagnostics

#ifdef _DEBUG
void CPlot::AssertValid() const
{
	CView::AssertValid();
}

void CPlot::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}
#endif //_DEBUG


// CPlot message handlers

void CPlot::OnDrawZoom()
{
	//Open the Zoom Dialog
	CZoom *zoomdialog;

	zoomdialog = new CZoom(&zoom,this);

	zoomdialog->Create(IDD_ZOOM,this);
}


void CPlot::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags) 
{
	short control;
	CPoint point = GetScrollPosition( );
	POINT pt;

	control = (GetKeyState(VK_CONTROL)&(-128));

      switch (nChar)
         {
         case VK_UP:
            //if (control) Cm_UP();
            //else {
				if (point.y!=0) {
					pt.x = point.x;
					pt.y = point.y - SCROLLCHANGE;
					ScrollToPosition(pt);
				}
            //}
            break;
         case VK_DOWN:
            //if (control) Cm_DOWN();
         	//else {
            	//if (point.y!=0) {
					pt.x = point.x;
					pt.y = point.y + SCROLLCHANGE;
					ScrollToPosition(pt);
				//}
            //}
            break;
         case VK_LEFT:
			 if (control) {
				if ((zoom - 10) >0) {
					zoom = zoom - 10;
					Invalidate();
				}

			 }
            else {
            	if (point.x!=0) {
					pt.x = point.x - SCROLLCHANGE;
					pt.y = point.y ;
					ScrollToPosition(pt);
				}
            }
            break;
         case VK_RIGHT:
			 if (control) {

				 zoom = zoom + 10;
				Invalidate();

			 }
            else {
            	//if (point.y!=0) {
					pt.x = point.x + SCROLLCHANGE;
					pt.y = point.y ;
					ScrollToPosition(pt);
				//}
            }
            break;
         case VK_HOME:
            //if (point.y!=0) {
					pt.x = 0;
					pt.y = 0;
					ScrollToPosition(pt);
				//}
            break;
         case VK_END:
            if (point.y!=0) {
					pt.x = point.x;
					pt.y = point.y - SCROLLCHANGE;
					ScrollToPosition(pt);
				}
            break;
         case VK_PRIOR:
            if (point.y!=0) {
					pt.x = point.x;
					pt.y = point.y + 10*SCROLLCHANGE;
					ScrollToPosition(pt);
				}
            break;
         case VK_NEXT:
            if (point.y!=0) {
					pt.x = point.x;
					pt.y = point.y - 10*SCROLLCHANGE;
					ScrollToPosition(pt);
				}
            break;

         //case '2':
         	//if (Scroller) Scroller->VScroll(SB_LINEDOWN,0);

         }
	
	CScrollView::OnKeyDown(nChar, nRepCnt, nFlags);
}


afx_msg void CPlot::OnLButtonDown( UINT nFlags, CPoint point ) {
	CString string;
	CPoint scroll;
	int i,j;
	char number[20];
	//check the location:
	int x,y;
	scroll = GetDeviceScrollPosition();
	x = point.x+scroll.x;
	y = point.y+scroll.y;
	
	i = (int) (floor) (((((float) y/(((float) zoom)/100.))-(float) topmargin-(float) extlinewidth))/((float)dotsize)) + plotdoc->ystart;
	j = (int) (floor) (((((float) x/(((float) zoom)/100.))-(float) leftmargin-(float) extlinewidth))/((float)dotsize)) + plotdoc->xstart;

	if (i>0&&i<=plotdoc->ct.GetSequenceLength()&&j>0&&j<=plotdoc->ct.GetSequenceLength()) {
		string = "BP ";
		itoa(i,number,10);
		string+=number;
		string+="(";
		string+=tobase(plotdoc->ct.numseq[i]);
		string+=")";
		string+="-";
		itoa(j,number,10);
		string+=number;
		string+="(";
		string+=tobase(plotdoc->ct.numseq[j]);
		string+=")";
		string+="  ;";
		string += plotdoc->message;
		string+="= ";
		if (plotdoc->arrayvalues[j][i]>INFINITE_ENERGY/10) {
			strcpy(number,plotdoc->outside);	
		}
		else gcvt(plotdoc->arrayvalues[j][i],6,number);
		string+=number;

		plotdoc->app->pMainFrame->m_wndStatusBar.SetWindowText(string);
	}
	
	

	CScrollView::OnLButtonDown(nFlags,point);
}

void CPlot::OnDrawColors()
{
	// Open A Dialog to Change the Number of Colors in Plot
	CColors *colordialog;

	colordialog = new CColors(plotdoc,this);

	colordialog->Create(IDD_COLORS,this);
}



void CPlot::OnDrawShowlegend()
{
	//Toggle the state of the "Show Legend Menu Item"
	CMenu* menu = plotdoc->app->pMainFrame->GetMenu();
	
	plotdoc->showlegend = !plotdoc->showlegend;
	if (plotdoc->showlegend) menu->CheckMenuItem(ID_DRAW_SHOWLEGEND,MF_CHECKED);
	else menu->CheckMenuItem(ID_DRAW_SHOWLEGEND,MF_UNCHECKED); 

	Invalidate();
}

void CPlot::OnDrawPlotrange()
{
	// Change the Range Displayed

	CPlotRange *plotrange;

	plotrange = new CPlotRange(plotdoc,this);

	plotrange->Create(IDD_PLOTRANGE,this);
}

void CPlot::OnFilePrint() 
{
	CView::OnFilePrint();
	
}

void CPlot::OnFilePrintPreview() 
{
	CView::OnFilePrintPreview();
	
}
BOOL CPlot::OnPreparePrinting(CPrintInfo* pInfo)
{
	// default preparation
	return DoPreparePrinting(pInfo);
}

void CPlot::OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo)
{
	//Scale to one page:
	pInfo->SetMaxPage(1);
}

void CPlot::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add cleanup after printing
}

void CPlot::OnPrint( CDC *pDC, CPrintInfo *pInfo ) {
	int pageheight, pagewidth;
	//float scalev,scaleh,scale;

	//Scale the image to fit one page:
	pageheight = pDC->GetDeviceCaps(VERTRES);
	pagewidth = pDC->GetDeviceCaps(HORZRES);
	//scalev = float(pageheight)/float(PageSize->cx+15);
	//scaleh = float(pagewidth)/float(PageSize->cy+15+(plotdoc->colors)*FontSize/2);

	//if(scalev<scaleh) scale = scalev;
	//else scale = scaleh;
	
	CSize  oldVExt, oldWExt;
	int prevMode;
	CSize sizeTotal;
	
	CRect rect;
	CPoint topleft,bottomright;
		
	prevMode = pDC->SetMapMode(MM_ISOTROPIC);
	
	

	oldWExt = pDC->SetWindowExt(WindowSize->cx,WindowSize->cy);
	oldVExt = pDC->SetViewportExt(pageheight,pagewidth);
	
	Draw(pDC);
	
	
	
	pDC->SetWindowExt(oldWExt);
	pDC->SetViewportExt(oldVExt);
	pDC->SetMapMode(prevMode);


}


void CPlot::OnOutputplot()
{
	//User has chosen to write the current plot to disk

	//Get the output filename:
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".ct","",OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"Dot Plot Files (*.dp)|*.dp|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {
		ofstream out;
		out.open(filedialog->GetPathName().GetBuffer(0));

		out << plotdoc->ct.GetSequenceLength()<<"\n";
		out << "i\tj\tkcal/mol\n";
		for (int i=plotdoc->xstart;i<=plotdoc->xstop;i++) {
			for (int j=max(i+1,plotdoc->ystart);j<=plotdoc->ystop;j++) {
				
				if (plotdoc->arrayvalues[j][i]<=plotdoc->colorrange[plotdoc->colors]&&plotdoc->arrayvalues[j][i]>=plotdoc->colorrange[0]) {
					out << i << "\t" << j << "\t" << plotdoc->arrayvalues[j][i] << "\n"; 
							
					
				}
				
			}

		}
				
		out.close();
	}
	delete filedialog;


	



}
