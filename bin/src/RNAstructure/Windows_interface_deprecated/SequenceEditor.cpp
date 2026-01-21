// SequenceEditor.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "SequenceEditor.h"

#include "RNAstructure.h"
#include <fstream>
#include <direct.h>
#include <mmsystem.h>

using namespace std;


#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#define REMAP 1101

/////////////////////////////////////////////////////////////////////////////
// CSequenceEditor




//This is the backend functions, set up to run in a seperate thread
UINT ReadProc( LPVOID pParam ){    
	CReadObject* pObject = (CReadObject*) pParam;
	CSequenceEditor* pEditor = (CSequenceEditor*) pObject->parent; 
	int i,j;
	char a[maxfil],c[maxfil],g[maxfil],u[maxfil],x[maxfil],t[maxfil];


	strcpy(a,pObject->datapath);
	strcat(a,"\\");
	strcat(a,"a.wav");

	strcpy(c,pObject->datapath);
	strcat(c,"\\");
	strcat(c,"c.wav");

	strcpy(g,pObject->datapath);
	strcat(g,"\\");
	strcat(g,"g.wav");

	strcpy(u,pObject->datapath);
	strcat(u,"\\");
	strcat(u,"u.wav");

	strcpy(t,pObject->datapath);
	strcat(t,"\\");
	strcat(t,"t.wav");

	strcpy(x,pObject->datapath);
	strcat(x,"\\");
	strcat(x,"x.wav");
	
	j = 0;
	for (i=pObject->start;i<=pObject->stop;i++) {
		if (pObject->cancel) break;
		
		if (pObject->texttoread[i]=='a'||pObject->texttoread[i]=='A') {
			
			
			::sndPlaySound(a, SND_SYNC);
			j++;
			Sleep(120);
			if (j%3==0) Sleep(400);
		}
		else if (pObject->texttoread[i]=='c'||pObject->texttoread[i]=='C') {
			
			
			::sndPlaySound(c, SND_SYNC);
			j++;
			Sleep(120);
			if (j%3==0) Sleep(400);
		}
		else if (pObject->texttoread[i]=='g'||pObject->texttoread[i]=='G') {
			
			//::PlaySound(a,NULL,SND_SYNC);
			::sndPlaySound(g, SND_SYNC);
			j++;
			Sleep(120);
			if (j%3==0) Sleep(400);
		}
		else if (pObject->texttoread[i]=='u'||pObject->texttoread[i]=='U') {
			
			//::PlaySound(a,NULL,SND_SYNC);
			::sndPlaySound(u, SND_SYNC);
			j++;
			Sleep(120);
			if (j%3==0) Sleep(400);
		}
		else if (pObject->texttoread[i]=='t'||pObject->texttoread[i]=='T') {
			
			//::PlaySound(a,NULL,SND_SYNC);
			::sndPlaySound(t, SND_SYNC);
			j++;
			Sleep(120);
			if (j%3==0) Sleep(400);
		}

		else if (pObject->texttoread[i]=='x'||pObject->texttoread[i]=='X') {
			
			//::PlaySound(a,NULL,SND_SYNC);
			::sndPlaySound(x, SND_SYNC);
			j++;
			Sleep(120);
			if (j%3==0) Sleep(400);
		}
		pObject->current = i;
		::SendMessage(pEditor->m_hWnd,ID_HIGHLIGHT,0,0);
		
		

	}
	

	
	::PostMessage(pEditor->m_hWnd,ID_READDONE,0,0);
	

	return 0;   // thread completed successfully
}





IMPLEMENT_DYNCREATE(CSequenceEditor, CFormView)

CSequenceEditor::CSequenceEditor()
	: CFormView(CSequenceEditor::IDD)
{
	//{{AFX_DATA_INIT(CSequenceEditor)
	m_Title = _T(" ");
	//}}AFX_DATA_INIT
	reading = false;
	readwhiletyping=false;
}

CSequenceEditor::~CSequenceEditor()
{
	delete SequenceWindow;
	delete CommentWindow;
	delete KpFont;

}


void CSequenceEditor::OnInitialUpdate() {
	CSequenceDoc *pDoc;
	CFont *pFont;
	
	pDoc = (CSequenceDoc*) GetDocument();
	SequenceWindow = new CSequenceEdit(&readwhiletyping,pDoc->datapath);
	CommentWindow = new CEdit();
	GetSeqDocument()->changed = false;
	
	ResizeParentToFit(FALSE);
	ResizeParentToFit();
	

	

	

	//Create the Sequence Editing Window
	CWnd* pWnd = GetDlgItem( IDC_SEQUENCE_POS );
	CRect rect;
	pWnd->GetWindowRect( &rect );
	ScreenToClient( &rect );

	
	SequenceWindow->Create(ES_MULTILINE|WS_VSCROLL   |WS_HSCROLL   |WS_VISIBLE|
		ES_WANTRETURN|WS_BORDER,rect,this,IDC_SEQUENCE_POS+REMAP);

	
	pWnd = GetDlgItem( IDC_COMMENT_POS );
	pWnd->GetWindowRect( &rect );
	ScreenToClient( &rect );

	
	CommentWindow->Create(ES_MULTILINE|WS_VSCROLL   |WS_HSCROLL   |WS_VISIBLE|
		ES_WANTRETURN|WS_BORDER,rect,this,IDC_COMMENT_POS+REMAP);	

	///////////////
	//pFont.CreateFont(12,0,0,0,FW_NORMAL,0,0,0,ANSI_CHARSET,OUT_TT_PRECIS,
	//	CLIP_TT_ALWAYS,PROOF_QUALITY,DEFAULT_PITCH|TMPF_TRUETYPE|FF_MODERN,"Courier New");
	//pFont.CreatePointFont(120,"Courier New");
	
	//pFont = GetFont();
	


	LOGFONT pLogFont;
	//pFont->GetLogFont(&pLogFont);


	pLogFont.lfHeight = 18;
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
	pLogFont.lfPitchAndFamily = FIXED_PITCH|FF_MODERN;
	strcpy(pLogFont.lfFaceName,"Courier New");

	KpFont = new CFont();
	KpFont->CreateFontIndirect(&pLogFont);


	//CommentWindow.SetFont(&pFont);
	SequenceWindow->SetFont(KpFont);
	//tFont(pFont);*/


	
	pFont = GetFont();
	CommentWindow->SetFont(pFont);

	//pFont.CreatePointFont(120,"Arial");
	//SequenceWindow.SetFont(&pFont);

	if (pDoc->filenamedefined) {
		//display the title, comment, and sequence

		m_Title=pDoc->title;
		SequenceWindow->SetWindowText(pDoc->sequence);
		CommentWindow->SetWindowText(pDoc->comment);
	}

	int lines,i,length;
	lines = SequenceWindow->GetLineCount();
	

	for (i=1;i<lines;i++) {
		length = SequenceWindow->LineLength(i);
		length++;

	}

	CView::OnInitialUpdate();

	
	
	UpdateData(FALSE);
}

void CSequenceEditor::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CSequenceEditor)
	DDX_Text(pDX, IDC_TITLE, m_Title);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CSequenceEditor, CFormView)
	ON_EN_CHANGE(IDC_SEQUENCE_POS+REMAP,OnChange)
	ON_EN_CHANGE(IDC_COMMENT_POS+REMAP,OnChange)
	ON_MESSAGE(ID_READDONE, OnReadDone)
	ON_MESSAGE(ID_HIGHLIGHT, OnHighLight)
	//{{AFX_MSG_MAP(CSequenceEditor)
	ON_COMMAND(ID_FILE_SAVESEQUENCE, OnFileSavesequence)
	ON_COMMAND(ID_FILE_SAVESEQUENCEAS, OnFileSavesequenceas)
	ON_EN_CHANGE(IDC_TITLE, OnChange)	
	ON_WM_DESTROY()
	ON_BN_CLICKED(IDC_FOLDASDNA, OnFoldasdna)
	ON_BN_CLICKED(IDC_FOLDASRNA, OnFoldasrna)
	ON_BN_CLICKED(IDC_READBACK, OnReadback)
	ON_COMMAND(ID_READ_READWHILETYPING, OnReadReadwhiletyping)
	ON_BN_CLICKED(IDC_FORMAT, Format)
	ON_COMMAND(ID_COPY, OnCopy)
	ON_COMMAND(ID_PASTE, OnPaste)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CSequenceEditor diagnostics

#ifdef _DEBUG
void CSequenceEditor::AssertValid() const
{
	CFormView::AssertValid();
}

void CSequenceEditor::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CSequenceEditor message handlers

void CSequenceEditor::OnFileSavesequence() 
{
	int i,j,length,lines;
	char *lineoftext,*lineoftext2;
	CSequenceDoc *pDoc;
	ofstream out;
	CString tempname;
	
		
	readobject.cancel = true;	


	//Save the sequence as a .seq file
	pDoc = (CSequenceDoc*) GetDocument();
	
	if (!pDoc->filenamedefined) {

		CFileDialog *filedialog;
		filedialog = new CFileDialog(FALSE,".seq",NULL,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
			"Sequence Files (*.seq)|*.seq|All Files|*.*||");
		filedialog->m_ofn.lpstrInitialDir=pDoc->startpath;
		if (filedialog->DoModal()==IDOK) {

			
			tempname=filedialog->GetPathName().GetBuffer(0);
			pDoc->allocatefilename(tempname.GetBuffer(10),tempname.GetLength()+1);
			//_getcwd(pDoc->startpath,_MAX_PATH);

			//now store the path in Startpath so that the program can start here next time:
			//_getcwd(GetFoldDocument()->startpath,_MAX_PATH);
			CString path;
			path = filedialog->GetPathName();
			i = path.GetLength();
			while(i>=0){
			
				if (path[i]=='\\') break;
				i--;
			}
			if (i>_MAX_PATH) i = _MAX_PATH;
			strncpy(pDoc->startpath,path.GetBuffer(1),i);
			*(pDoc->startpath + i) ='\0';


			delete filedialog;

		}
		else {
			
			delete filedialog;
			return;
		}


	}

	out.open(pDoc->filename);

	

	lines = CommentWindow->GetLineCount();
	

	
	length = CommentWindow->GetWindowTextLength();
	lineoftext = new char [length+1];
	CommentWindow->GetWindowText(lineoftext,length+1);
	
	lineoftext2 = new char [length+1+2*lines];

	strcpy(lineoftext2,";");
	j = 0;
	for (i=0;i<length+1;i++) {
		if (lineoftext[i]=='\r') {
			strncat(lineoftext2,lineoftext+j,i-j);
			strcat(lineoftext2,"\n;");
			
			if(lineoftext[i+1]=='\n') {
				i++;
					
			}
			j=i+1;

		}
		else if (lineoftext[i]=='\n') {
			strncat(lineoftext2,lineoftext+j,i-j);
			strcat(lineoftext2,"\n;");
			j=i+1;

		}

	}
	if (j!=i+1) {
		strcat(lineoftext2,lineoftext+j);
		strcat(lineoftext2,"\n");
	}
	else lineoftext2[strlen(lineoftext2)-2]='\0';

	out<<lineoftext2;
	//out << "\n";


	delete[] lineoftext;
	delete[] lineoftext2;
	UpdateData(TRUE);
	
	out << m_Title.GetBuffer(10);
	out << "\n";

	/*lines = SequenceWindow->GetLineCount();
	

	for (i=0;i<lines;i++) {
		length = SequenceWindow->LineLength(i);
		if (length>0) {
			lineoftext = new char[length+1];
			SequenceWindow->GetLine(i,lineoftext,length);
			lineoftext[length]='\0';
	
			out <<lineoftext;
			out <<"\n";
		
			delete[] lineoftext;
		}
		else {
			out << "\n";
		}

	}*/

	lines = SequenceWindow->GetLineCount();
	

	
	length = SequenceWindow->GetWindowTextLength();
	lineoftext = new char [length+1];
	SequenceWindow->GetWindowText(lineoftext,length+1);
	
	lineoftext2 = new char [length+1+2*lines];

	strcpy(lineoftext2,"");
	j = 0;
	for (i=0;i<length+1;i++) {
		if (lineoftext[i]=='\r') {
			strncat(lineoftext2,lineoftext+j,i-j);
			strcat(lineoftext2,"\n");
			
			if(lineoftext[i+1]=='\n') {
				i++;
					
			}
			j=i+1;

		}
		else if (lineoftext[i]=='\n') {
			strncat(lineoftext2,lineoftext+j,i-j);
			strcat(lineoftext2,"\n");
			j=i+1;

		}

	}
	if (j!=i+1) {
		strcat(lineoftext2,lineoftext+j);
		
	}
	out <<lineoftext2;
	out << "1\n";

	delete[] lineoftext;
	delete[] lineoftext2;

	out.close();

	GetSeqDocument()->changed = false;
	return;

	
}

CSequenceDoc* CSequenceEditor::GetSeqDocument() {


	return (CSequenceDoc*) GetDocument();
	
}

void CSequenceEditor::OnFileSavesequenceas() 
{
	CSequenceDoc *pDoc;
	CString tempname;
	readobject.cancel = true;
	pDoc = GetSeqDocument(); 

	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".seq",NULL,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"Sequence Files (*.seq)|*.seq|All Files|*.*||");
	filedialog->m_ofn.lpstrInitialDir=pDoc->startpath;
	if (filedialog->DoModal()==IDOK) {

			
		tempname=filedialog->GetPathName().GetBuffer(0);
		pDoc->allocatefilename(tempname.GetBuffer(10),tempname.GetLength()+1);
		//_getcwd(pDoc->startpath,_MAX_PATH);

		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(GetFoldDocument()->startpath,_MAX_PATH);
		CString path;
		path = filedialog->GetPathName();
		int i = path.GetLength();
		while(i>=0){
			
			if (path[i]=='\\') break;
			i--;
		}
		if (i>_MAX_PATH) i = _MAX_PATH;
		strncpy(pDoc->startpath,path.GetBuffer(1),i);
		*(pDoc->startpath + i) ='\0';

		delete filedialog;

		OnFileSavesequence();
	}
	else {
			
		delete filedialog;
			
	}
	
}

void CSequenceEditor::Format()
{
	readobject.cancel = true;
	SequenceWindow->format();
/*	char **lines;
	//CString *cstringlines;
	int linecount,i,j,count;
	CString output;
	int foo;

	readobject.cancel = true;

	linecount = SequenceWindow->GetLineCount();
	lines = new char *[linecount];
	//cstringlines = new CString [linecount];
	

	for (i=1;i<linecount;i++) {
		foo = ((CEdit*) SequenceWindow)->LineLength(i)+1;
		lines[i] = new char[SequenceWindow->LineLength(i)+1];
		SequenceWindow->GetLine(i,lines[i]);
		lines[i][SequenceWindow->LineLength(i)]='\0';
		foo=foo;

		//SequenceWindow->GetLine(i,cstringlines[i].GetBuffer(1));
		//lines[i] = new char [cstringlines[i].GetLength()];
		//strcpy(lines[i],cstringlines[i].GetBuffer(1));


	}

	output = "";

	count = 0;
	for (i=0;i<linecount;i++) {
		for (j=0;j<=((int)strlen(lines[i]));j++) {

			if(lines[i][j]=='A'||lines[i][j]=='C'||lines[i][j]=='G'||
				lines[i][j]=='U'||lines[i][j]=='T'||lines[i][j]=='a'||
				lines[i][j]=='c'||lines[i][j]=='g'||lines[i][j]=='u'||
				lines[i][j]=='t'||lines[i][j]=='x'||lines[i][j]=='X'||
				lines[i][j]=='I'||lines[i][j]=='i') {

				output+=lines[i][j];
				count++;

				if(count%50==0) {
					output+="\r\n";
				}
				else if (count%5==0) {

					output+=" ";
				}

			}


		}

	}


	SequenceWindow->SetWindowText(output);
	//delete[] cstringlines;
	for (i=0;i<linecount;i++) delete[] lines[i];
	delete[] lines;*/



}



void CSequenceEditor::OnChange() 
{
	//This function is called anytime the text in a field is changed:
	GetSeqDocument()->changed = true;
	
}






void CSequenceEditor::OnDestroy() 
{
	readobject.cancel = true;
	if (GetSeqDocument()->changed) {

		if (!ShouldSave()) return;

	}
	
	
	CFormView::OnDestroy();
	
	
	
}

void CSequenceEditor::OnFoldasdna() 
{
	// User Clicked the Fold as DNA button
	readobject.cancel = true;
	Fold(false);
	
}

void CSequenceEditor::OnFoldasrna() 
{
	// User Clicked the Fold as RNA button
	readobject.cancel = true;
	Fold(true);
	
}


void CSequenceEditor::Fold(bool isRNA) {

	CSequenceDoc *pDoc;
	CString tempname;

	readobject.cancel = true;

	pDoc = (CSequenceDoc*) GetDocument();


	if (GetSeqDocument()->changed) {
		if (!ShouldSave()) return;
		GetSeqDocument()->changed = false;


	}

	if (isRNA) 
		((CRNAstructureApp*) pDoc->pMainFrame)->FoldRNA(pDoc->filename);
	else 
		((CRNAstructureApp*) pDoc->pMainFrame)->FoldDNA(pDoc->filename);
	

	GetParentFrame()->SendMessage(WM_CLOSE);

}

bool CSequenceEditor::ShouldSave() {
	int i;
	

	i = AfxMessageBox( "The sequence has been edited.\nSave changes?", 
			MB_YESNOCANCEL   |MB_ICONQUESTION   );

	if (i==IDCANCEL) return false;
	else if (i==IDYES) OnFileSavesequence();


	return true;

}
void CSequenceEditor::OnReadback() 
{

	
	readobject.cancel = true;	
	SequenceWindow->SetFocus();
if (!reading) {
	reading = true;	
	
		
	//Start a thread for reading back
	
	//First get the window text:
	readobject.texttoread = new char [SequenceWindow->GetWindowTextLength()+1];
	SequenceWindow->GetWindowText(readobject.texttoread,SequenceWindow->GetWindowTextLength()+1);
	readobject.cancel = false;

	SequenceWindow->GetSel( readobject.start, readobject.stop);
	if (readobject.start==NULL) readobject.start = 0;
	if (readobject.stop==NULL) readobject.stop = SequenceWindow->GetWindowTextLength();
	else if (readobject.start==readobject.stop) {
		readobject.start = 0;
		SequenceWindow->GetWindowTextLength();
	}
	else readobject.stop--;
	
	readobject.datapath = GetSeqDocument()->datapath;
	readobject.parent = this;

	CButton* button = (CButton*) GetDlgItem( IDC_READBACK );
	button->SetWindowText("Cancel Read");

	SequenceWindow->SetFocus();

	AfxBeginThread(ReadProc,&readobject);

}

	
}

LRESULT CSequenceEditor::OnReadDone(WPARAM wParam, LPARAM lParam) {

	delete[] readobject.texttoread;
	reading = false;
	CButton* button = (CButton*) GetDlgItem( IDC_READBACK );
	button->SetWindowText("Read Sequence");
	return 0;
}

LRESULT CSequenceEditor::OnHighLight(WPARAM i, LPARAM j) {
	SequenceWindow->SetSel( readobject.current+1, readobject.current+2);
	return 0;
}

void CSequenceEditor::OnReadReadwhiletyping() 
{
	readwhiletyping=(!readwhiletyping);

	CMenu* menu = GetSeqDocument()->pFrame->GetMenu( );
	
	
	if (readwhiletyping) menu->CheckMenuItem(ID_READ_READWHILETYPING,MF_CHECKED);
	else menu->CheckMenuItem(ID_READ_READWHILETYPING,MF_UNCHECKED); 


	
	
}

void CSequenceEditor::OnCopy() 
{
	SequenceWindow->Copy();
	
}

void CSequenceEditor::OnPaste() 
{
	SequenceWindow->Paste();
	
}
