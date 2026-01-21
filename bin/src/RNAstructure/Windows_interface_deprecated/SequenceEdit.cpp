// SequenceEdit.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "SequenceEdit.h"
#include <mmsystem.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CSequenceEdit

IMPLEMENT_DYNCREATE(CSequenceEdit, CEdit)

CSequenceEdit::CSequenceEdit(bool *Readwhiletyping,char *datapath)
{

	read = Readwhiletyping;

	


	strcpy(a,datapath);
	strcat(a,"\\");
	strcat(a,"a.wav");

	strcpy(c,datapath);
	strcat(c,"\\");
	strcat(c,"c.wav");

	strcpy(g,datapath);
	strcat(g,"\\");
	strcat(g,"g.wav");

	strcpy(u,datapath);
	strcat(u,"\\");
	strcat(u,"u.wav");

	strcpy(t,datapath);
	strcat(t,"\\");
	strcat(t,"t.wav");

	strcpy(x,datapath);
	strcat(x,"\\");
	strcat(x,"x.wav");
}

CSequenceEdit::CSequenceEdit()
{


}

CSequenceEdit::~CSequenceEdit()
{
}


BEGIN_MESSAGE_MAP(CSequenceEdit, CEdit)
	//{{AFX_MSG_MAP(CSequenceEdit)
	ON_WM_KEYDOWN()
	ON_WM_CHAR()
	ON_WM_SYSCHAR()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CSequenceEdit drawing



/////////////////////////////////////////////////////////////////////////////
// CSequenceEdit diagnostics

#ifdef _DEBUG
void CSequenceEdit::AssertValid() const
{
	CEdit::AssertValid();
}

void CSequenceEdit::Dump(CDumpContext& dc) const
{
	CEdit::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CSequenceEdit message handlers

void CSequenceEdit::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags) 
{
	
	short control;
	control = (GetKeyState(VK_CONTROL)&(-128));
	if ((nChar=='V'||nChar=='v')) {
		//The user pressed cntrl-V
		if (control) Paste();

	}

	
	CEdit::OnKeyDown(nChar, nRepCnt, nFlags);
}

void CSequenceEdit::OnSysChar(UINT nChar, UINT nRepCnt, UINT nFlags) {

	CWnd::OnSysChar(nChar, nRepCnt, nFlags);

}

void CSequenceEdit::OnChar(UINT nChar, UINT nRepCnt, UINT nFlags) 
{

	short control;
	control = (GetKeyState(VK_CONTROL)&(-128));
	

	if ((nChar=='V'||nChar=='v')) {
		//The user pressed cntrl-V
		if (control) Paste();

	}
	// Only allow certain characters
	if (nChar == 65||nChar == 67||nChar == 71||nChar == 85||
		nChar == 84||nChar == 97||nChar == 99||nChar == 103||
		nChar == 117||nChar == 116||nChar == 88||nChar == 120||
		nChar == 73||nChar == 65||nChar == 9||
		nChar == 32||nChar == 8||nChar == 13) {
		CEdit::OnChar(nChar, nRepCnt, nFlags);

		if (*read) {
			switch (nChar) {
				case (65): 
					::sndPlaySound(a, SND_SYNC);
					break;
				case (84): 
					::sndPlaySound(t, SND_SYNC);
					break;
				case (67): 
					::sndPlaySound(c, SND_SYNC);
					break;
				case (71): 
					::sndPlaySound(g, SND_SYNC);
					break;
				case (85): 
					::sndPlaySound(u, SND_SYNC);
					break;
				case (97): 
					::sndPlaySound(a, SND_SYNC);
					break;
				case (99): 
					::sndPlaySound(c, SND_SYNC);
					break;
				case (103): 
					::sndPlaySound(g, SND_SYNC);
					break;
				case (117): 
					::sndPlaySound(u, SND_SYNC);
					break;
				case (116): 
					::sndPlaySound(t, SND_SYNC);
					break;
				case (88): 
					::sndPlaySound(x, SND_SYNC);
					break;
				case (120): 
					::sndPlaySound(x, SND_SYNC);
					break;
				
				


			}
		}

	}

	
}

void CSequenceEdit::format() {
	//char **lines;
	
	//int linecount,i,j,count;
	CString output;
	CString input;
	//int foo;
	//UpdateData(true);
	int count,i;
	

	//linecount = GetLineCount();
	//lines = new char *[linecount];
	
	
	GetWindowText(input);
	//for (i=0;i<linecount;i++) {
	//	foo = LineLength(i+1);
	//	lines[i] = new char[LineLength(i+1)+5];
	//	GetLine(i,lines[i]);
	//	lines[i][LineLength(i+1)+1]='\0';
	//	foo =foo;

		//SequenceWindow->GetLine(i,cstringlines[i].GetBuffer(1));
		//lines[i] = new char [cstringlines[i].GetLength()];
		//strcpy(lines[i],cstringlines[i].GetBuffer(1));


	//}



	output = "";

	count = 0;
	for (i=0;i<input.GetLength();i++) {
		/*for (j=0;j<=((int)strlen(lines[i]));j++) {

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


		}*/

		if (input[i]=='A'||input[i]=='C'||input[i]=='G'||
				input[i]=='U'||input[i]=='T'||input[i]=='a'||
				input[i]=='c'||input[i]=='g'||input[i]=='u'||
				input[i]=='t'||input[i]=='x'||input[i]=='X'||
				input[i]=='I'||input[i]=='i') {

			output+=input[i];
			count++;

			if(count%50==0) {
				output+="\r\n";
			}
			else if (count%5==0) {

				output+=" ";
			}

		}

	}


	SetWindowText(output);
	//delete[] cstringlines;
	//for (i=0;i<linecount;i++) {
	//	delete[] lines[i];
	//}
	//delete[] lines;
}