// SequenceDoc2.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "SequenceDoc.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CSequenceDoc2

IMPLEMENT_DYNCREATE(CSequenceDoc, CDocument)

CSequenceDoc::CSequenceDoc()
{
}

CSequenceDoc::CSequenceDoc(CMainFrame* Frame, char *Startpath, char *Datapath, CWinApp *MainFrame, char *Filename)
{

	startpath = Startpath;
	datapath = Datapath;
	if (Filename!=NULL) {
		filename = new char [strlen(Filename)+1];
		strcpy(filename,Filename);
		filenamedefined=true;

		SetTitle(filename);

		//Now open the file as well
		check = OpenSequence(filename);
	}
	else filenamedefined=false;

	pMainFrame = MainFrame;
	pFrame = Frame;

}

BOOL CSequenceDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;
	return TRUE;
}

CSequenceDoc::~CSequenceDoc()
{
	if (filenamedefined) delete[] filename;

}

void CSequenceDoc::allocatefilename(char *name,int length) {


	filename = new char [length+1];
	strcpy(filename,name);
	filenamedefined = true;
	SetTitle(filename);

}




BEGIN_MESSAGE_MAP(CSequenceDoc, CDocument)
	//{{AFX_MSG_MAP(CSequenceDoc2)
		// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CSequenceDoc2 diagnostics

#ifdef _DEBUG
void CSequenceDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CSequenceDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CSequenceDoc2 serialization

void CSequenceDoc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}
}

/////////////////////////////////////////////////////////////////////////////
// CSequenceDoc2 commands

int CSequenceDoc::OpenSequence(char *filename) {

//Open a sequence:
	//This routine recognizes Genbank and mfold style sequences

	ifstream in;
	char line[maxsequencelinelength],line2[maxsequencelinelength];
	bool genbank,plain,standard;
	int i,length;
	FILE *se;





	//etstablish the sequence type:
	genbank = false;
	plain = false;
	standard = false;
	in.open(filename);

	in >> line;
	if(line[0]!=';') {
		while(!in.eof()) {

			in >> line;
			if (!strcmp(line,"ORIGIN")) {
				//it's a Genbank file
				genbank = true;
				break;
			}
			else {
				for (i=0;i<((int)strlen(line));i++) {
					if (!(line[i]=='A'||line[i]=='C'||line[i]=='G'||line[i]=='U'||
						line[i]=='T'||line[i]=='a'||line[i]=='c'||line[i]=='g'||
						line[i]=='u'||line[i]=='t'||line[i]=='x'||line[i]=='X'||
						line[i]=='i'||line[i]=='I'||line[i]==' '||line[i]=='\n')) {

						//something other than a sequence indicator was found
						plain = false;

					}

				}

			}


		}



	}
	else {
		//Its a standard sequence file
		standard = true;

	}

	in.close();

	if (genbank) {
		
		
		//This time read the sequence
		genbank = false;
		in.open(filename);
		in >> line;
		while(!in.eof()&&!genbank) {
			if (!strcmp("DEFINITION",line)) {

				comment = ""; 
				in >> line;
				while (!strcmp(line,"ACCESSION")) {

					comment += line;

				}
				

			}
			else if (!strcmp("LOCUS",line)) {

				title="";
				in >> line;
				title+=line;

			}
			else if (!strcmp(line,"ORIGIN")) {
				//got to the sequence portion:
				genbank = true;//use this to break from while statement
				length = 0;
				sequence="";
				in >> line;
				while(strcmp(line,"//")) {
					for (i=0;i<((int)strlen(line));i++) {
						if (line[i]=='A'||line[i]=='C'||line[i]=='G'||line[i]=='U'||
						line[i]=='T'||line[i]=='x'||line[i]=='X') {
							
							length++;
							sequence += line[i];
							if (length%50==0) sequence+="\n\r";
							else if (length%5==0) sequence+="\t";



						}
						else if (line[i]=='a') {
							length++;
							sequence += "A";
							if (length%50==0) sequence+="\n\r";
							else if (length%5==0) sequence+="\t";

						}
						else if (line[i]=='c') {
							length++;
							sequence += "C";
							if (length%50==0) sequence+="\n\r";
							else if (length%5==0) sequence+="\t";

						}
						else if (line[i]=='g') {

							length++;
							sequence += "G";
							if (length%50==0) sequence+="\n\r";
							else if (length%5==0) sequence+="\t";
						}
						else if (line[i]=='u') {
							length++;
							sequence += "U";
							if (length%50==0) sequence+="\n\r";
							else if (length%5==0) sequence+="\t";

						}
						else if (line[i]=='t') {
							length++;
							sequence += "T";
							if (length%50==0) sequence+="\n\r";
							else if (length%5==0) sequence+="\t";

						}

						else if (line[i]=='N'||line[i]=='n') {
							length++;
							sequence += "X";
							if (length%50==0) sequence+="\n\r";
							else if (length%5==0) sequence+="\t";

						}


					}

					in >> line;
				}


				
			}
			in >> line;


		}


		in.close();
	}
	else if (plain) {
		
		length = 0;
		in.open(filename);
		comment ="";
		title="";
		sequence="";

		in >> line;
		while (in.eof()) {
			for (i=0;i<((int)strlen(line));i++) {
				if (line[i]=='A'||line[i]=='C'||line[i]=='G'||line[i]=='U'||
						line[i]=='T'||line[i]=='a'||line[i]=='c'||line[i]=='g'||
						line[i]=='u'||line[i]=='t'||line[i]=='x'||line[i]=='X'||
						line[i]=='i'||line[i]=='I') {
							
							length++;
							sequence+=line[i];
							if (length%50==0) sequence+="\r\n";
							else if (length%5==0) sequence+="\t";
				}

			}
			in >> line;

		}


		in.close();


	}
	else if (standard) {

		//read the sequence file to get the number of nucleotides
		se=fopen(filename,"r");
		title ="";
		comment="";
		sequence="";

		fgets(line,maxsequencelinelength,se);
		while (line[0]==';') {
			strcpy(line2,line+1);
			line2[strlen(line2)-1] = '\0';
			
			comment += line2;
			

			fgets(line,maxsequencelinelength,se);
			if(line[0]==';') comment+="\r\n";
		}
		line[strlen(line)-1]='\0';
		title+=line;
		
		while (standard){
			
			fgets(line,maxsequencelinelength,se);

			if (line[strlen(line)-2]=='1') {
				standard = false;
				line[strlen(line)-2]='\0';
				
				
				sequence += line;



			}
			else if (line[strlen(line)-1]=='1') {
				standard = false;
				line[strlen(line)-1]='\0';
				
				
				sequence += line;



			}
			else {
				line[strlen(line)-1]='\0';
				sequence+=line;
				sequence+="\r\n";

			}
		}
		
		


	}
	else return 0;

	return 1;


}
