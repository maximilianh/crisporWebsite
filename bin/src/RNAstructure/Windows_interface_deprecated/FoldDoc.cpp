// FoldDoc.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "FoldDoc.h"
#include "../src/defines.h"
#include "../src/algorithm.h"
#include "globals.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CFoldDoc

IMPLEMENT_DYNCREATE(CFoldDoc, CCTDoc)


//isRNA = true is RNA parameters should be loaded
//		= false for DNA parameters
CFoldDoc::CFoldDoc(bool isRNA, char *datapath, CWinApp *pMainFrame1, char* startpath1, bool fold, bool *Savefile, bool PF)
{	
	char loop[maxfil],stackf[maxfil], tstackh[maxfil], tstacki[maxfil],
		tloop[maxfil], miscloop[maxfil], danglef[maxfil], int22[maxfil],
		int21[maxfil],coax[maxfil], tstackcoax[maxfil],coaxstack[maxfil],
		tstack[maxfil], tstackm[maxfil],triloop[maxfil],int11[maxfil],hexaloop[maxfil],
		tstacki23[maxfil],tstacki1n[maxfil];

		ISRNA = isRNA;
		Datapath = datapath;

		
	 
		getdat(loop, stackf, tstackh, tstacki,
			tloop, miscloop, danglef, int22,
			int21,coax, tstackcoax,
			coaxstack, tstack, tstackm, triloop,
			int11, hexaloop, tstacki23, tstacki1n, datapath, isRNA);

		if (opendat(loop, stackf, tstackh, tstacki,
			tloop, miscloop, danglef, int22,
			int21,coax, tstackcoax,
			coaxstack, tstack, tstackm, triloop,
			int11,hexaloop,tstacki23, tstacki1n, &data)) OK = true;
		else OK = false;

		pMainFrame = pMainFrame1;

		startpath = startpath1;

		if (!PF) {
			if (isRNA&&fold) {
				SetTitle("RNA Fold");
			}
			else if (fold) SetTitle("DNA Fold");
			else if (isRNA) SetTitle("RNA Efn2");
			else SetTitle("DNA Efn2");
		}
		else {
			if (isRNA&&fold) {
				SetTitle("RNA Partition Function");
			}
			else if (fold) SetTitle("DNA Partition Function");

		}
	
		filenamedefined = false;

		savefile = Savefile;

		T = (float) 310.15;

		maximuminternalloopsize = 30;



}

int CFoldDoc::newtemp() {
	//The temperature of folding was changed from 37 degrees C, alter the data files:
	char loop[maxfil],stackf[maxfil], tstackh[maxfil], tstacki[maxfil],
		tloop[maxfil], miscloop[maxfil], danglef[maxfil], int22[maxfil],
		int21[maxfil],coax[maxfil], tstackcoax[maxfil],coaxstack[maxfil],
		tstack[maxfil], tstackm[maxfil],triloop[maxfil],int11[maxfil],hexaloop[maxfil],
		tstacki23[maxfil],tstacki1n[maxfil];
	datatable *enthalpy;

	getenthalpydat(loop, stackf, tstackh, tstacki,
			tloop, miscloop, danglef, int22,
			int21,coax, tstackcoax,
			coaxstack, tstack, tstackm, triloop,
			int11, hexaloop, tstacki23, tstacki1n, Datapath, ISRNA);

	enthalpy = new datatable();

	if (opendat(loop, stackf, tstackh, tstacki,
			tloop, miscloop, danglef, int22,
			int21,coax, tstackcoax,
			coaxstack, tstack, tstackm, triloop,
			int11,hexaloop,tstacki23, tstacki1n, enthalpy)==0) {
	
				delete enthalpy;
				return 0;

	}
	dG_T(T,data,*enthalpy,data);

	delete enthalpy;

	return 1;

}

void CFoldDoc::namefile(char *Filename) {
	int i;

	filenamedefined = true;
	
	i = strlen(Filename);

	filename = new char [i+2];
	strcpy (filename,Filename);


}

BOOL CFoldDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;
	return TRUE;
}

CFoldDoc::~CFoldDoc()
{
	if(filenamedefined) delete[] filename;


}


BEGIN_MESSAGE_MAP(CFoldDoc, CCTDoc)
	//{{AFX_MSG_MAP(CFoldDoc)
		// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CFoldDoc diagnostics

#ifdef _DEBUG
void CFoldDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CFoldDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CFoldDoc serialization

void CFoldDoc::Serialize(CArchive& ar)
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
// CFoldDoc commands
