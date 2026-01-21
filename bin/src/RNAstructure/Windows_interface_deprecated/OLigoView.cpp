// OLigoView.cpp : implementation file
//


#include "stdafx.h"
//#include <string>
#include "RNAstructure.h"
#include "OligoView.h"
#include "godialog.h"
#include "oligodialog.h"
#include "../src/algorithm.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#define PAIRED_COLOR 255
#define OVERALL_G_COLOR 16711680
#define DUPLEX_G_COLOR 7633158
#define BREAKTARGET_G_COLOR 7633158
#define INTRAMOLEC_G_COLOR 16711935
#define INTERMOLEC_G_COLOR 4194432
#define GREEN_COLOR RGB(0,255,0)

#define BUTTONWIDTH 50
#define BUTTONHEIGHT 25



/////////////////////////////////////////////////////////////////////////////
// COLigoView

IMPLEMENT_DYNCREATE(COligoView, CView)

COligoView::COligoView()
{

	LOGFONT pLogFont;

	current = 1;
	

	fontsize = 24;

	//set the initial sequence size
	sequencesize = 20;


	sequence = new CFont();
	pLogFont.lfHeight = sequencesize;
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
	sequence->CreateFontIndirect(&pLogFont);


	symbol = new CFont();
	pLogFont.lfHeight = fontsize;
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
	strcpy(pLogFont.lfFaceName,"Symbol");
	symbol->CreateFontIndirect(&pLogFont);


	subscript = new CFont();
	pLogFont.lfHeight = fontsize/2;
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
	subscript->CreateFontIndirect(&pLogFont);

	font = new CFont();
	pLogFont.lfHeight = fontsize;
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
	font->CreateFontIndirect(&pLogFont);

	poverallBrush.CreateSolidBrush(OVERALL_G_COLOR);
	predbrush.CreateSolidBrush(PAIRED_COLOR);
	potherBrush = new CBrush();
	
	
	
}

COligoView::~COligoView()
{

	delete GetDocument()->oligoobject;
	delete symbol;
	delete subscript;
	delete font;

	delete go;
	delete left;
	delete right;
	delete tenleft;
	delete tenright;

	delete potherBrush;
}


BEGIN_MESSAGE_MAP(COligoView, CView)
	ON_BN_CLICKED(ID_GO, OnGo)
	ON_BN_CLICKED(ID_LEFT, OnLeft)
	ON_BN_CLICKED(ID_TENLEFT, OnTenLeft)
	ON_BN_CLICKED(ID_RIGHT, OnRight)
	ON_BN_CLICKED(ID_TENRIGHT, OnTenRight)
	ON_WM_KEYDOWN()
	ON_WM_LBUTTONDOWN()
	//{{AFX_MSG_MAP(COLigoView)
	ON_COMMAND(ID_GRAPH_FREEENERGY_OVERALL, OnGraphFreeenergyOverall)
	ON_COMMAND(ID_GRAPH_FREEENERGY_OVERALLANDDUPLEX, OnGraphFreeenergyOverallandduplex)
	ON_COMMAND(ID_GRAPH_FREEENERGY_BROKENTARGETSTRUCTURE, OnGraphFreeenergyBrokentargetstructure)
	ON_COMMAND(ID_GRAPH_FREEENERGY_DUPLEX, OnGraphFreeenergyDuplex)
	ON_COMMAND(ID_GRAPH_FREEENERGY_OLIGOMERBIMOLECULAR, OnGraphFreeenergyOligomerbimolecular)
	ON_COMMAND(ID_GRAPH_FREEENERGY_OLIGOMERUNIMOLECULAR, OnGraphFreeenergyOligomerunimolecular)
	ON_WM_LBUTTONDBLCLK()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// COLigoView drawing

void COligoView::OnInitialUpdate()
{
	
	RECT rect,rectwin;
	GetClientRect(&rect);
	

	 

	//put the buttons on the screen
	go = new CButton();
	left  = new CButton();
	right = new CButton();
	tenleft = new CButton();
	tenright = new CButton();

	rectwin.left = rect.right/2;
	rectwin.right = rectwin.left+BUTTONWIDTH;
	rectwin.top = fontsize*4;
	rectwin.bottom = fontsize*4+BUTTONHEIGHT;
	go->Create("Go...", BS_PUSHBUTTON|WS_CHILD|WS_VISIBLE|WS_TABSTOP  ,rectwin, this, ID_GO);

	rectwin.left = (long)(rect.right/2+(1.5)*BUTTONWIDTH);
	rectwin.right = rectwin.left+BUTTONWIDTH;
	right->Create(">", BS_PUSHBUTTON|WS_CHILD|WS_VISIBLE|WS_TABSTOP  ,rectwin, this, ID_RIGHT);

	rectwin.left = rect.right/2+(3)*BUTTONWIDTH;
	rectwin.right = rectwin.left+BUTTONWIDTH;
	tenright->Create(">>", BS_PUSHBUTTON|WS_CHILD|WS_VISIBLE|WS_TABSTOP  ,rectwin, this, ID_TENRIGHT);

	rectwin.left = (long)(rect.right/2-(1.5)*BUTTONWIDTH);
	rectwin.right = rectwin.left+BUTTONWIDTH;
	left->Create("<", BS_PUSHBUTTON|WS_CHILD|WS_VISIBLE|WS_TABSTOP  ,rectwin, this, ID_LEFT);

	rectwin.left = rect.right/2-(3)*BUTTONWIDTH;
	rectwin.right = rectwin.left+BUTTONWIDTH;
	tenleft->Create("<<", BS_PUSHBUTTON|WS_CHILD|WS_VISIBLE|WS_TABSTOP  ,rectwin, this, ID_TENLEFT);
	 
	
	CView::OnInitialUpdate();
	OnGraphFreeenergyOverallandduplex(); 
	
	//CSize sizeTotal;
	// TODO: calculate the total size of this view
	//sizeTotal.cx = sizeTotal.cy = 100;
	//SetScrollSizes(MM_TEXT, sizeTotal);
}

int COligoView::TextOut(CString text,CDC* pDC, int x, int y, bool draw) {
	CSize size;

	if (draw) {
		pDC->TextOut(x,y,text);
	}
	size = pDC->GetTextExtent(text);
	return size.cx;

}

int COligoView::DGOut(CDC* pDC, int x, int y, bool draw) {
	
	int totalwidth,width;
	//std::string temperature; 
	char temperature[10];
	

	
	CFont *oldFont = pDC->SelectObject(symbol);
	
	
	if (draw) {
		width = TextOut("D",pDC,x,y,true);
		pDC->SelectObject(oldFont);
		totalwidth=width;
		width = TextOut("G",pDC,x+width,y,true);
		totalwidth+=width;
		TextOut("° = ",pDC,x+totalwidth,y,true);

		

		
		pDC->SelectObject(subscript);

		
		//temperature = (GetDocument()->oligoobject->T-273.15);
		sprintf(temperature,"%.1f",GetDocument()->oligoobject->T-273.15);

		TextOut(temperature,pDC,x+totalwidth,y+fontsize/2,true);
		pDC->SelectObject(oldFont);
		
		return 0;	

	}
	else {
		//determine the size and return it
		
		totalwidth=TextOut("D", pDC, x,y,false);
		width = totalwidth;
		pDC->SelectObject(oldFont);
		totalwidth = TextOut("G° = ",pDC,x,y,false);
		width+=totalwidth;
		
		return width;

	}

	
	

}

void COligoView::OnDraw(CDC* pDC)
{
	COligoDoc* pDoc = GetDocument();
	
	CSize size;
	int width,line2,line3,midline,line4,line5,column3,i,
		columnwidth,graphtop,graphbottom,graphzero;
	RECT rect;
	char number[10];
	double cp;
	double dg;
	
	CPen *pen;


	GetClientRect(&rect);
	midline = rect.right/2;

	column2 = 2*rect.right/5;
	column3 = 4*rect.right/5;

	
	CFont *oldFont = pDC->SelectObject(font);

	
	pDC->TextOut(0,0,"Oligo #");
	size = pDC->GetTextExtent("O");
	line2 = size.cy+fontsize/5;
	line3 = 2*line2;
	line4 = 3*line2;
	line5 = 4*line2;
	oligoline = 6*line2;
	sequenceline = 6*line2+sequencesize;

	itoa(current,number,10);
	pDC->TextOut(size.cy,line2,number);

	if (GetDocument()->oligoobject->isdna) pDC->TextOut(0,line3,"DNA Oligo");
	else pDC->TextOut(0,line3,"RNA Oligo");

	if ((GetDocument()->oligoobject->c)<=1e-12) {
		cp = (GetDocument()->oligoobject->c)/1e-12;
		_gcvt(cp,8,number);
		pDC->TextOut(0,line4,number);
		size = pDC->GetTextExtent(number);
		width = size.cx;
		size = pDC->GetTextExtent(" ");
		width+=size.cx;
		pDC->TextOut(width,line4,"pM");
	}
	else if ((GetDocument()->oligoobject->c)<=1e-9) {
		cp = (GetDocument()->oligoobject->c)/1e-9;
		_gcvt(cp,8,number);
		pDC->TextOut(0,line4,number);
		size = pDC->GetTextExtent(number);
		width = size.cx;
		size = pDC->GetTextExtent(" ");
		width+=size.cx;
		pDC->TextOut(width,line4,"nM");
	}
	else if ((GetDocument()->oligoobject->c)<=1e-6) {
		cp = (GetDocument()->oligoobject->c)/1e-6;
		_gcvt(cp,8,number);
		pDC->TextOut(0,line4,number);
		size = pDC->GetTextExtent(number);
		width = size.cx;
		size = pDC->GetTextExtent(" ");
		width+=size.cx;
		pDC->TextOut(width,line4,"uM");
	}
	else {
		
		cp = (GetDocument()->oligoobject->c)/1e-3;
		_gcvt(cp,8,number);
		pDC->TextOut(0,line4,number);
		size = pDC->GetTextExtent(number);
		width = size.cx;
		size = pDC->GetTextExtent(" ");
		width+=size.cx;
		pDC->TextOut(width,line4,"mM");

	}


	//put the DG info on the screen:
	
	if (overallandduplex||overall) pDC->SetTextColor(OVERALL_G_COLOR);
	//overall:
	dg = ((float) (GetDocument()->oligoobject->table[current][0]))/10.0;
	gcvt(dg,5,number);
	pDC->TextOut(column2,0,number);
	width = DGOut(pDC,column2,0,false);
	DGOut(pDC,column2-width,0,true);

	width+=TextOut("Overall ",pDC,column2-width,0,false);
	TextOut("Overall",pDC,column2-width,0,true);
	if (overallandduplex||overall) pDC->SetTextColor(0);

	if (overallandduplex||overall) pDC->SetTextColor(DUPLEX_G_COLOR);
	
	//duplex
	dg = ((float) ( GetDocument()->oligoobject->table[current][1]))/10.0;
	gcvt(dg,5,number);
	pDC->TextOut(column2,line2,number);
	width = DGOut(pDC,column2,line2,false);
	DGOut(pDC,column2-width,line2,true);

	width+=TextOut("Duplex ",pDC,column2-width,line2,false);
	TextOut("Duplex",pDC,column2-width,line2,true);

	if (overallandduplex||overall) pDC->SetTextColor(0);
	

	//Tm
	dg = ((float) (GetDocument()->oligoobject->table[current][5]))/10.0;
	gcvt(dg,5,number);
	pDC->TextOut(column2,line3,number);
	width = TextOut("Tm = ",pDC,column2,line3,false);
	width = TextOut("Tm =",pDC,column2-width,line3,true);

	if (target) pDC->SetTextColor(BREAKTARGET_G_COLOR );
	//target
	dg = ((float) (GetDocument()->oligoobject->table[current][2]))/10.0;
	gcvt(dg,5,number);
	pDC->TextOut(column3,0,number);
	width = DGOut(pDC,column3,0,false);
	DGOut(pDC,column3-width,0,true);

	width+=TextOut("Break targ. ",pDC,column3-width,0,false);
	TextOut("Break targ.",pDC,column3-width,0,true);

	if (target) pDC->SetTextColor(0);


	if (intramolecular) pDC->SetTextColor(INTRAMOLEC_G_COLOR);
	//intramolecular
	dg = ((float) (GetDocument()->oligoobject->table[current][3]))/10.0;
	gcvt(dg,5,number);
	pDC->TextOut(column3,line2,number);
	width = DGOut(pDC,column3,line2,false);
	DGOut(pDC,column3-width,line2,true);

	width+=TextOut("oligo-self ",pDC,column3-width,line2,false);
	TextOut("oligo-self",pDC,column3-width,line2,true);
	if (intramolecular) pDC->SetTextColor(0);

	if (intermolecular) pDC->SetTextColor(INTERMOLEC_G_COLOR);
	//intermolecular
	dg = ((float) (GetDocument()->oligoobject->table[current][4]))/10.0;
	gcvt(dg,5,number);
	pDC->TextOut(column3,line3,number);
	width = DGOut(pDC,column3,line3,false);
	DGOut(pDC,column3-width,line3,true);

	width+=TextOut("oligo-oligo ",pDC,column3-width,line3,false);
	TextOut("oligo-oligo",pDC,column3-width,line3,true);
	if (intermolecular) pDC->SetTextColor(0);
	

	

	//Now Draw the sequence
	
	//put the oligo with the left edge at column2 on oligoline
	

	pDC->SelectObject(sequence);
	width = TextOut("'",pDC,column2,oligoline,false);
	TextOut("'",pDC,column2-width,oligoline,true);
	width = (int)(1.5*(double)width);
	TextOut("3",pDC,column2-width,oligoline,true);


	//If siRNA design is being performed, color the oligo green for acceptible and red for not
	if (GetDocument()->oligoobject->siRNA) {
	
		if (GetDocument()->oligoobject->mask[current]) pDC->SetTextColor(GREEN_COLOR);
		else pDC->SetTextColor(PAIRED_COLOR);

	}

	width = column2;
	for (i=current;i<current+GetDocument()->oligoobject->length;i++) {
		
		switch(GetDocument()->oligoobject->ct.numseq[i]) {
			case (0) :
				width+=TextOut("X",pDC,width,oligoline,true);
				break;

			case(1):

				if (GetDocument()->oligoobject->isdna) width+=TextOut("T",pDC,width,oligoline,true);
				else width+=TextOut("U",pDC,width,oligoline,true);
				break;

			case(2):
				width+=TextOut("G",pDC,width,oligoline,true);
				break;

			case(3):
				width+=TextOut("C",pDC,width,oligoline,true);
				break;


			case(4):
				width+=TextOut("A",pDC,width,oligoline,true);
				break;


			case(5):
				width+=TextOut("X",pDC,width,oligoline,true);
				break;

		}
	}
	if (GetDocument()->oligoobject->siRNA) {
		pDC->SetTextColor(0);

	}

	i=TextOut("5",pDC,width,oligoline,false);
	TextOut("'",pDC,(int)(width+0.5*i),oligoline,true);
	TextOut("5",pDC,width,oligoline,true);

	//Now drop the target sequence:
	i =current;
	width = column2;
	seqtextwidth = TextOut("X",pDC,width,sequenceline,false);
	while (width+TextOut("G",pDC,width,sequenceline,false)<rect.right&&
		i<=GetDocument()->oligoobject->ct.GetSequenceLength()) {

		if (GetDocument()->oligoobject->ct.GetPair(i,1)>0) pDC->SetTextColor(PAIRED_COLOR);
		switch(GetDocument()->oligoobject->ct.numseq[i]) {
			case (0) :
				
				width+=TextOut("X",pDC,width,sequenceline,true);
				
				break;

			case(1):
				width+=TextOut("A",pDC,width,sequenceline,true);
				
				break;

			case(2):
				width+=TextOut("C",pDC,width,sequenceline,true);
				break;

			case(3):
				width+=TextOut("G",pDC,width,sequenceline,true);
				break;


			case(4):
				
				if (GetDocument()->oligoobject->istargetdna) width+=TextOut("T",pDC,width,oligoline,true);
				else width+=TextOut("U",pDC,width,sequenceline,true);
				break;


			case(5):
				width+=TextOut("I",pDC,width,sequenceline,true);
				break;

		}
		if (GetDocument()->oligoobject->ct.GetPair(i,1)>0) pDC->SetTextColor(0);
		i++;

	}
	
	i = current-1;
	width = column2-TextOut("G",pDC,width,sequenceline,false);

	while (width>0&&i>0) {

		if (GetDocument()->oligoobject->ct.GetPair(i,1)>0) pDC->SetTextColor(PAIRED_COLOR);
		switch(GetDocument()->oligoobject->ct.numseq[i]) {
			case (0) :
				width-=TextOut("X",pDC,width,sequenceline,true);
				break;

			case(1):
				width-=TextOut("A",pDC,width,sequenceline,true);
				
				break;

			case(2):
				width-=TextOut("C",pDC,width,sequenceline,true);
				break;

			case(3):
				width-=TextOut("G",pDC,width,sequenceline,true);
				break;


			case(4):
				
				if (GetDocument()->oligoobject->istargetdna) width-=TextOut("T",pDC,width,oligoline,true);
				else width-=TextOut("U",pDC,width,sequenceline,true);
				break;


			case(5):
				width-=TextOut("I",pDC,width,sequenceline,true);
				break;

		}
		if (GetDocument()->oligoobject->ct.GetPair(i,1)>0) pDC->SetTextColor(0);
		i--;

	}


	//Now Draw the bars for the graph
	
	graphtop = sequenceline+2*fontsize;
	graphbottom = rect.bottom - fontsize;
	graphzero = graphtop+(int)((float)(graphbottom-graphtop))*(((float)(graphmax))/((float) (graphmax-graphmin)));

	i =current;
	columnwidth =TextOut("G",pDC,width,sequenceline,false); 
	width = column2;
	
	pen = new CPen;
	COLORREF black = 0;
	pen->CreatePen(PS_SOLID,1,black);
	pDC->SelectObject(pen);

	while (width+columnwidth<rect.right&&
		i<=GetDocument()->oligoobject->ct.GetSequenceLength()-GetDocument()->oligoobject->length+1) {
	
		
		if (duplex||overallandduplex){
			DrawBar(pDC, i, 1, width, columnwidth, graphzero, graphtop, graphbottom, potherBrush);

		}
		
		if (overall||overallandduplex) {
			DrawBar(pDC, i, 0, width, columnwidth, graphzero, graphtop, graphbottom, &poverallBrush);

		}

		if (target) {

			DrawBar(pDC, i, 2, width, columnwidth, graphzero, graphtop, graphbottom, potherBrush);
		}
		
		if (intramolecular) {

			DrawBar(pDC, i, 3, width, columnwidth, graphzero, graphtop, graphbottom, potherBrush);
		}

		if (intermolecular) {

			DrawBar(pDC, i, 4, width, columnwidth, graphzero, graphtop, graphbottom, potherBrush);
		}
		
		/*duplex=false;
		target=false;
		duplex=false;
		intramolecular=false;
		intermolecular=true;*/
		
		
		width+=columnwidth;
		i++;
	}

	int end = width; 
	
	
	

	i = current-1;
	width = column2-columnwidth;

	//Draw a scale and a line for zero
	int largest;
	

	largest = TextOut("0",pDC,width,sequenceline,false);
	
	dg = ((float) (graphmax)/10.0);
	gcvt(dg,5,number);
	largest = max(largest, TextOut(number,pDC,width,sequenceline,false));
	dg = ((float) (graphmin)/10.0);
	gcvt(dg,5,number);
	largest = max(largest, TextOut(number,pDC,width,sequenceline,false));



	while (width>largest&&i>0) {

	
		
		if (duplex||overallandduplex){
			DrawBar(pDC, i, 1, width, columnwidth, graphzero, graphtop, graphbottom, potherBrush);

		}
		if (overall||overallandduplex) {
			DrawBar(pDC, i, 0, width, columnwidth, graphzero, graphtop, graphbottom, &poverallBrush);
		}
	
		
		if (target) {

			DrawBar(pDC, i, 2, width, columnwidth, graphzero, graphtop, graphbottom, potherBrush);
		}
		
		if (intramolecular) {

			DrawBar(pDC, i, 3, width, columnwidth, graphzero, graphtop, graphbottom, potherBrush);
		}

		if (intermolecular) {

			DrawBar(pDC, i, 4, width, columnwidth, graphzero, graphtop, graphbottom, potherBrush);
		}
		

		width-=columnwidth;
		i--;

	}

		//Draw the current bar in red:

		if (duplex){
			DrawBar(pDC, current, 1, column2, columnwidth, graphzero, graphtop, graphbottom, &predbrush);

		}
		if (overall||overallandduplex) {
			DrawBar(pDC, current, 0, column2, columnwidth, graphzero, graphtop, graphbottom, &predbrush);
		}
	
		
		if (target) {

			DrawBar(pDC, current, 2, column2, columnwidth, graphzero, graphtop, graphbottom, &predbrush);
		}
		
		if (intramolecular) {

			DrawBar(pDC, current, 3, column2, columnwidth, graphzero, graphtop, graphbottom, &predbrush);
		}

		if (intermolecular) {

			DrawBar(pDC, current, 4, column2, columnwidth, graphzero, graphtop, graphbottom, &predbrush);
		}

	//Draw a scale and a line for zero
	
	width= width+columnwidth;

	pDC->MoveTo(width,graphtop);
	pDC->LineTo(width,graphbottom);

	pDC->MoveTo(width-2,graphzero);
	pDC->LineTo(end,graphzero);
	
	pDC->MoveTo(width-2,graphtop);
	pDC->LineTo(end,graphtop);

	pDC->MoveTo(width-2,graphbottom);
	pDC->LineTo(end,graphbottom);

	TextOut("0",pDC,width-TextOut("0",pDC,width,sequenceline,false)-2,graphzero-fontsize/2,true);
	dg = ((float) (graphmax)/10.0);
	gcvt(dg,5,number);
	TextOut(number,pDC,width-TextOut(number,pDC,width,sequenceline,false)-2,graphtop-fontsize/2,true);
	dg = ((float) (graphmin)/10.0);
	gcvt(dg,5,number);
	TextOut(number,pDC,width-TextOut(number,pDC,width,sequenceline,false)-2,graphbottom-fontsize/2,true);

	delete pen;


	//restore the DC
	
	pDC->SelectObject(oldFont);
	
}

/////////////////////////////////////////////////////////////////////////////
// COLigoView diagnostics

#ifdef _DEBUG
void COligoView::AssertValid() const
{
	CView::AssertValid();
}

void COligoView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// COLigoView message handlers




COligoDoc *COligoView::GetDocument() {

	return ((COligoDoc*) CView::GetDocument());

}

void COligoView::OnLeft() {

	current--;
	if (current<1) current =1;
	InvalidateRect(NULL);
	GetDocument()->Frame->UpdateWindow();

}

void COligoView::OnRight() {

	current++;
	if (current>GetDocument()->oligoobject->ct.GetSequenceLength()-GetDocument()->oligoobject->length+1) 
		current=GetDocument()->oligoobject->ct.GetSequenceLength()-GetDocument()->oligoobject->length+1;
	InvalidateRect(NULL);
	UpdateWindow();

}

void COligoView::OnGo() {
	CGoDialog *Go;

	Go = new CGoDialog(&current, GetDocument()->oligoobject, this);
	Go->Create(IDD_GODIALOG,this);

}

void COligoView::OnTenLeft() {

	current-=10;
	if (current<1) current =1;
	InvalidateRect(NULL);
	UpdateWindow();

}

void COligoView::OnTenRight(){

	current+=10;
	if (current>GetDocument()->oligoobject->ct.GetSequenceLength()-GetDocument()->oligoobject->length+1) 
		current=GetDocument()->oligoobject->ct.GetSequenceLength()-GetDocument()->oligoobject->length+1;
	InvalidateRect(NULL);
	UpdateWindow();

}


void COligoView::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags) 
{
	short control;
	

	control = (GetKeyState(VK_CONTROL)&(-128));

      switch (nChar)
         {
         
            
         case VK_LEFT:
         	if (control) OnTenLeft();
            else {
            	OnLeft();
            }
            break;
         case VK_RIGHT:
         	if (control) OnTenRight();
            else {
            	OnRight();
            }
            break;
         case VK_HOME:
            current = 1;
			InvalidateRect(NULL);
			UpdateWindow();


            break;
         case VK_END:
            current=GetDocument()->oligoobject->ct.GetSequenceLength()-GetDocument()->oligoobject->length+1;
			InvalidateRect(NULL);
			UpdateWindow();
            break;
         

         }
	
	
}


void COligoView::OnGraphFreeenergyOverall() 
{
	int i;
	
	//Put check on menu item
	overallandduplex=false;
	overall=true;
	duplex=false;
	target=false;
	duplex=false;
	intramolecular=false;
	intermolecular=false;
	
	placecheck();
	
	graphmin = 0;
	graphmax = 0;
	//Get some variables ready for drawing the graph
	for (i=1;i<=GetDocument()->oligoobject->ct.GetSequenceLength()-GetDocument()->oligoobject->length+1;i++) {
		if (GetDocument()->oligoobject->table[i][0]<graphmin) graphmin = GetDocument()->oligoobject->table[i][0];
		else if (GetDocument()->oligoobject->table[i][0]>graphmax) graphmax = GetDocument()->oligoobject->table[i][0];
	}
	
	InvalidateRect(NULL);
	GetDocument()->Frame->UpdateWindow();	
	
}

void COligoView::OnGraphFreeenergyOverallandduplex() 
{
	int i;
	
	//Put check on menu item
	overallandduplex=true;
	overall=false;
	duplex=false;
	target=false;
	duplex=false;
	intramolecular=false;
	intermolecular=false;
	
	placecheck();
	
	graphmin = 0;
	graphmax = 0;
	//Get some variables ready for drawing the graph
	for (i=1;i<=GetDocument()->oligoobject->ct.GetSequenceLength()-GetDocument()->oligoobject->length+1;i++) {
		if (GetDocument()->oligoobject->table[i][0]<graphmin) graphmin = GetDocument()->oligoobject->table[i][0];
		else if (GetDocument()->oligoobject->table[i][0]>graphmax) graphmax = GetDocument()->oligoobject->table[i][0];
		if (GetDocument()->oligoobject->table[i][1]<graphmin) graphmin = GetDocument()->oligoobject->table[i][1];
		
		else if (GetDocument()->oligoobject->table[i][1]>graphmax) graphmax = GetDocument()->oligoobject->table[i][1];
	}
	delete potherBrush;
	potherBrush = new CBrush();
	potherBrush->CreateSolidBrush(DUPLEX_G_COLOR);

	InvalidateRect(NULL);
	GetDocument()->Frame->UpdateWindow();

	
}

void COligoView::placecheck() {
	CMenu* menu = GetDocument()->pmainframe->GetMenu( );
	
	
	menu->CheckMenuItem(ID_GRAPH_FREEENERGY_OVERALLANDDUPLEX,MF_UNCHECKED);
	menu->CheckMenuItem(ID_GRAPH_FREEENERGY_OVERALL,MF_UNCHECKED);
	menu->CheckMenuItem(ID_GRAPH_FREEENERGY_DUPLEX,MF_UNCHECKED);
	menu->CheckMenuItem(ID_GRAPH_FREEENERGY_BROKENTARGETSTRUCTURE,MF_UNCHECKED);
	menu->CheckMenuItem(ID_GRAPH_FREEENERGY_DUPLEX,MF_UNCHECKED);
	menu->CheckMenuItem(ID_GRAPH_FREEENERGY_OLIGOMERUNIMOLECULAR,MF_UNCHECKED);
	menu->CheckMenuItem(ID_GRAPH_FREEENERGY_OLIGOMERBIMOLECULAR,MF_UNCHECKED);

	if (overallandduplex) menu->CheckMenuItem(ID_GRAPH_FREEENERGY_OVERALLANDDUPLEX,MF_CHECKED);
	else if (overall) menu->CheckMenuItem(ID_GRAPH_FREEENERGY_OVERALL,MF_CHECKED);
	else if (duplex) menu->CheckMenuItem(ID_GRAPH_FREEENERGY_DUPLEX,MF_CHECKED);
	else if (target) menu->CheckMenuItem(ID_GRAPH_FREEENERGY_BROKENTARGETSTRUCTURE,MF_CHECKED);
	else if (duplex) menu->CheckMenuItem(ID_GRAPH_FREEENERGY_DUPLEX,MF_CHECKED);
	else if (intramolecular) menu->CheckMenuItem(ID_GRAPH_FREEENERGY_OLIGOMERUNIMOLECULAR,MF_CHECKED);
	else if (intermolecular) menu->CheckMenuItem(ID_GRAPH_FREEENERGY_OLIGOMERBIMOLECULAR,MF_CHECKED);


}

void COligoView::OnGraphFreeenergyBrokentargetstructure() 
{
	int i;
	
	//Put check on menu item
	overallandduplex=false;
	overall=false;
	duplex=false;
	target=true;
	duplex=false;
	intramolecular=false;
	intermolecular=false;
	
	placecheck();
	
	graphmin = 0;
	graphmax = 0;
	//Get some variables ready for drawing the graph
	for (i=1;i<=GetDocument()->oligoobject->ct.GetSequenceLength()-GetDocument()->oligoobject->length+1;i++) {
		if (GetDocument()->oligoobject->table[i][2]<graphmin) graphmin = GetDocument()->oligoobject->table[i][2];
		else if (GetDocument()->oligoobject->table[i][2]>graphmax) graphmax = GetDocument()->oligoobject->table[i][2];
	}
	delete potherBrush;
	potherBrush = new CBrush();
	potherBrush->CreateSolidBrush(BREAKTARGET_G_COLOR);

	InvalidateRect(NULL);
	GetDocument()->Frame->UpdateWindow();


}

void COligoView::OnGraphFreeenergyDuplex() 
{
	int i;
	
	//Put check on menu item
	overallandduplex=false;
	overall=false;
	duplex=true;
	target=false;
	intramolecular=false;
	intermolecular=false;
	
	placecheck();
	
	graphmin = 0;
	graphmax = 0;
	//Get some variables ready for drawing the graph
	for (i=1;i<=GetDocument()->oligoobject->ct.GetSequenceLength()-GetDocument()->oligoobject->length+1;i++) {
		if (GetDocument()->oligoobject->table[i][1]<graphmin) graphmin = GetDocument()->oligoobject->table[i][1];
		else if (GetDocument()->oligoobject->table[i][1]>graphmax) graphmax = GetDocument()->oligoobject->table[i][1];
	}

	delete potherBrush;
	potherBrush = new CBrush();
	potherBrush->CreateSolidBrush(DUPLEX_G_COLOR);

	InvalidateRect(NULL);
	GetDocument()->Frame->UpdateWindow();

	
}

void COligoView::OnGraphFreeenergyOligomerbimolecular() 
{
	int i;
	
	//Put check on menu item
	overallandduplex=false;
	overall=false;
	duplex=false;
	target=false;
	duplex=false;
	intramolecular=false;
	intermolecular=true;
	
	placecheck();
	
	graphmin = 0;
	graphmax = 0;
	//Get some variables ready for drawing the graph
	for (i=1;i<=GetDocument()->oligoobject->ct.GetSequenceLength()-GetDocument()->oligoobject->length+1;i++) {
		if (GetDocument()->oligoobject->table[i][4]<graphmin) graphmin = GetDocument()->oligoobject->table[i][4];
		else if (GetDocument()->oligoobject->table[i][4]>graphmax) graphmax = GetDocument()->oligoobject->table[i][4];
	}
	delete potherBrush;
	potherBrush = new CBrush();
	potherBrush->CreateSolidBrush(INTERMOLEC_G_COLOR);
	InvalidateRect(NULL);
	GetDocument()->Frame->UpdateWindow();
}

void COligoView::OnGraphFreeenergyOligomerunimolecular() 
{
	int i;
	
	//Put check on menu item
	overallandduplex=false;
	overall=false;
	duplex=false;
	target=false;
	duplex=false;
	intramolecular=true;
	intermolecular=false;
	
	placecheck();
	
	graphmin = 0;
	graphmax = 0;
	//Get some variables ready for drawing the graph
	for (i=1;i<=GetDocument()->oligoobject->ct.GetSequenceLength()-GetDocument()->oligoobject->length+1;i++) {
		if (GetDocument()->oligoobject->table[i][3]<graphmin) graphmin = GetDocument()->oligoobject->table[i][3];
		else if (GetDocument()->oligoobject->table[i][3]>graphmax) graphmax = GetDocument()->oligoobject->table[i][3];
	}
	delete potherBrush;
	potherBrush = new CBrush();
	potherBrush->CreateSolidBrush(INTRAMOLEC_G_COLOR);
	InvalidateRect(NULL);
	GetDocument()->Frame->UpdateWindow();	
}


int COligoView::findheight(int graphzero,int bar,int graphtop, int graphbottom) {


	return graphtop+(int)((float)(graphbottom-graphtop))*(((float)(graphmax - bar))/((float) (graphmax-graphmin)));

}

afx_msg void COligoView::OnLButtonDown( UINT nFlags, CPoint point ) {
	
	//check the location and change current:
	if (point.y>sequenceline) {

		
		current = current + (point.x-column2)/seqtextwidth;
		if (point.x<column2) current--;
		if (current<1) current = 1;
		else if (current>GetDocument()->oligoobject->ct.GetSequenceLength()-GetDocument()->oligoobject->length+1)
			current = GetDocument()->oligoobject->ct.GetSequenceLength()-GetDocument()->oligoobject->length+1;
		InvalidateRect(NULL);
		GetDocument()->Frame->UpdateWindow();


	}


	CView::OnLButtonDown(nFlags,point);
}


void COligoView::OnLButtonDblClk(UINT nFlags, CPoint point) 
{
	structure ct;
	bool uni;
	COligoDialog oligodialog;
	int i;
	char character[100];

	if (point.y>oligoline&&point.y<sequenceline) {
		
		//fold and draw the current oligo
		if (GetDocument()->oligoobject->table[current][3]<0&&GetDocument()->oligoobject->table[current][4]<0) {
			//ask the user whether to show the bi or uni molecular folding
			if (oligodialog.DoModal()==IDOK) {
				uni = true;
				
			}
			else {
				uni=false;
				
			}

		}
		else if (GetDocument()->oligoobject->table[current][4]<0) {
			
			uni = false;
		}
		else {
			uni=true;
			
		}

		int slen;
		if (uni) {
			
			slen = GetDocument()->oligoobject->length;
		}
		else slen = 2*GetDocument()->oligoobject->length+3;

		ct.allocate(slen);
		string name;
		name = "Oligowalk oligomer #";
		//strcpy(ct.ctlabel[1],);
		itoa(current,character,10);
		//strcat(ct.ctlabel[1],character);
		name+=character;
		name+="\n";
		//strcat(ct.ctlabel[1],"\n");
		ct.SetSequenceLabel(name);
		for (i=1;i<=GetDocument()->oligoobject->length;i++) {
			ct.hnumber[i]=i;
			ct.numseq[ct.GetSequenceLength()-i+1]=complement(current+i-1,&GetDocument()->oligoobject->ct);
			ct.nucs[ct.GetSequenceLength()-i+1] = (tobase(ct.numseq[ct.GetSequenceLength()-i+1]))[0]; 
			if (!uni) {
				ct.numseq[(ct.GetSequenceLength()-3)/2-i+1]=complement(current+i-1,&GetDocument()->oligoobject->ct); 	
				ct.nucs[(ct.GetSequenceLength()-3)/2-i+1] =(tobase(ct.numseq[(ct.GetSequenceLength()-3)/2-i+1]))[0];
				
			}
		}
		if (!uni) {
			ct.numseq[GetDocument()->oligoobject->length+1]=5;
			ct.numseq[GetDocument()->oligoobject->length+2]=5;
			ct.numseq[GetDocument()->oligoobject->length+3]=5;
			ct.intermolecular = true;
			ct.inter[0]=GetDocument()->oligoobject->length+1;
			ct.inter[1]=GetDocument()->oligoobject->length+2;
			ct.inter[2]=GetDocument()->oligoobject->length+3;
			ct.nucs[GetDocument()->oligoobject->length+1] ='I';
			ct.nucs[GetDocument()->oligoobject->length+2] ='I';
			ct.nucs[GetDocument()->oligoobject->length+3] ='I';

		}

		if (GetDocument()->oligoobject->isdna) {
			dynamic(&ct,GetDocument()->oligoobject->ddata,100,20,0);
			//efn2 (GetDocument()->oligoobject->ddata,&ct);
		}	
		else {
			dynamic(&ct,&GetDocument()->oligoobject->data,100,20,0);
			//efn2 (&GetDocument()->oligoobject->data,&ct);

		}
			
			
		//sortstructures (&ct);
		//ct.sort();
		ct.ctout("%%oligo_walk_oligomer.ct");
		((CRNAstructureApp*) GetDocument()->pParent)->Draw("%%oligo_walk_oligomer.ct");



	}
	
	CView::OnLButtonDblClk(nFlags, point);
}


void COligoView::DrawBar(CDC *pDC, int i, int array, int x, int width, int graphzero, int graphtop, int graphbottom, CBrush *pBrush) {
	int height;

	CRgn *pRgn;
	pRgn = new CRgn();
	height=findheight(graphzero,GetDocument()->oligoobject->table[i][array],graphtop,graphbottom);
	if (GetDocument()->oligoobject->table[i][array]>0) {
		pRgn->CreateRectRgn(x, height, 
		x+width,graphzero );
	}
	else {
		pRgn->CreateRectRgn(x,graphzero,  
			x+width ,height );

	}
			


	pDC->FillRgn( pRgn, pBrush );
	pDC->MoveTo(x,height);
	pDC->LineTo(x,graphzero);
	pDC->LineTo(x+width,graphzero);
	pDC->LineTo(x+width,height);
	pDC->LineTo(x,height);
	delete pRgn;

}
