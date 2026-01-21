// CTDoc.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "CTDoc.h"

#include "../src/algorithm.h"

#include <iostream>
using namespace std;

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CCTDoc

IMPLEMENT_DYNCREATE(CCTDoc, CDocument)

CCTDoc::CCTDoc()
{
}

BOOL CCTDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;
	return TRUE;
}

CCTDoc::~CCTDoc()
{
}


int CCTDoc::integer(char nuc) {

	if (nuc=='A'||nuc=='a') return 1;
	else if (nuc=='c'||nuc=='C') return 2;
	else if (nuc=='g'||nuc=='G') return 3;
	else if (nuc=='U'||nuc=='u'||nuc=='T'||nuc=='t') return 4;
	
	else if (nuc=='i'||nuc=='I') return 5;
	else return 0;
	


}


int CCTDoc::OpenSequence(char *filename) {
	//Open a sequence:
	//This routine recognizes Genbank and mfold style sequences

	ifstream in;
	char line[maxsequencelinelength];
	bool genbank,plain;
	int i,length;


	//etstablish the sequence type:
	genbank = false;
	plain = false;
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
		in.close();
		return ct.openseq(filename);

	}

	in.close();

	if (genbank) {
		//get the length of the sequence
		genbank = false;
		in.open(filename);
		in >> line;
		while(!in.eof()&&!genbank) {
			if (!strcmp(line,"ORIGIN")) {
				//got to the sequence portion:
				genbank = true;//use this to break from while statement
				length = 0;
				in >> line;
				while(strcmp(line,"//")) {
					for (i=0;i<((int)strlen(line));i++) {
						if (line[i]=='A'||line[i]=='C'||line[i]=='G'||line[i]=='U'||
						line[i]=='T'||line[i]=='a'||line[i]=='c'||line[i]=='g'||
						line[i]=='u'||line[i]=='t'||line[i]=='x'||line[i]=='X'||
						line[i]=='i'||line[i]=='I') {
							length++;	

						}


					}

					in >> line;
				}


				
			}


		}


		in.close();

		ct.allocate(length);
		//This time read the sequence
		genbank = false;
		in.open(filename);
		in >> line;
		while(!in.eof()&&!genbank) {
			if (!strcmp("DEFINITION",line)) {

				//strcpy(ct.ctlabel[1],"");
				in >> line;
				string label;
				label="";
				while (!strcmp(line,"ACCESSION")) {
					label+=line;
					label+=" ";
					//strcat(ct.ctlabel[1],line);
					//strcat(ct.ctlabel[1]," ");

				}
				label+="\n";
				ct.SetSequenceLabel(label);
				//strcat(ct.ctlabel[1],"\n");

			}
			else if (!strcmp(line,"ORIGIN")) {
				//got to the sequence portion:
				genbank = true;//use this to break from while statement
				length = 0;
				in >> line;
				while(strcmp(line,"//")) {
					for (i=0;i<((int)strlen(line));i++) {
						if (line[i]=='A'||line[i]=='C'||line[i]=='G'||line[i]=='U'||
						line[i]=='T'||line[i]=='a'||line[i]=='c'||line[i]=='g'||
						line[i]=='u'||line[i]=='t'||line[i]=='x'||line[i]=='X'||
						line[i]=='i'||line[i]=='I') {
							
							length++;
							ct.numseq[length] = integer(line[i]);
							
							ct.nucs[i] = line[i];



						}


					}

					in >> line;
				}


				
			}


		}


		in.close();
	}
	else if (plain) {
		in.open(filename); 
		length = 0;
		in >> line;
		while (in.eof()) {
			for (i=0;i<((int)strlen(line));i++) {
				if (line[i]=='A'||line[i]=='C'||line[i]=='G'||line[i]=='U'||
						line[i]=='T'||line[i]=='a'||line[i]=='c'||line[i]=='g'||
						line[i]=='u'||line[i]=='t'||line[i]=='x'||line[i]=='X'||
						line[i]=='i'||line[i]=='I') {
							
							length++;
				}

			}
			in >> line;

		}


		in.close();

		ct.allocate(length);
		length = 0;
		in.open(filename);

		in >> line;
		while (in.eof()) {
			for (i=0;i<((int)strlen(line));i++) {
				if (line[i]=='A'||line[i]=='C'||line[i]=='G'||line[i]=='U'||
						line[i]=='T'||line[i]=='a'||line[i]=='c'||line[i]=='g'||
						line[i]=='u'||line[i]=='t'||line[i]=='x'||line[i]=='X'||
						line[i]=='i'||line[i]=='I') {
							
							length++;
							ct.numseq[length] = integer(line[i]);
							
							ct.nucs[i] = line[i];
				}

			}
			in >> line;

		}


		in.close();


	}
	else return 0;

	return 1;

}


BEGIN_MESSAGE_MAP(CCTDoc, CDocument)
	//{{AFX_MSG_MAP(CCTDoc)
		// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CCTDoc diagnostics

#ifdef _DEBUG
void CCTDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CCTDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CCTDoc serialization

void CCTDoc::Serialize(CArchive& ar)
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
// CCTDoc commands
