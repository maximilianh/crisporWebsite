// FoldDoc.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "DynDoc.h"
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

IMPLEMENT_DYNCREATE(CDynDoc, CCTDoc)


//isRNA = true is RNA parameters should be loaded
//		= false for DNA parameters
CDynDoc::CDynDoc(bool isRNA, char *datapath, CWinApp *pMainFrame1, char* startpath1, bool Singleinsert)
{	
	char loop[maxfil],stackf[maxfil], tstackh[maxfil], tstacki[maxfil],
		tloop[maxfil], miscloop[maxfil], danglef[maxfil], int22[maxfil],
		int21[maxfil],coax[maxfil], tstackcoax[maxfil],coaxstack[maxfil],
		tstack[maxfil], tstackm[maxfil],triloop[maxfil],int11[maxfil],hexaloop[maxfil],
		tstacki23[maxfil],tstacki1n[maxfil];

		ISRNA = isRNA;
	 
		getdat(loop, stackf, tstackh, tstacki,
			tloop, miscloop, danglef, int22,
			int21,coax, tstackcoax,
			coaxstack, tstack, tstackm, triloop,
			int11, hexaloop, tstacki23, tstacki1n, datapath, isRNA);


		//dynalign requires an aletrnative filename for miscloop.dat
		//strcpy(miscloop,datapath);
		//strcat(miscloop,"\\");
		//if (isRNA) strcat(miscloop,"dynalignmiscloop.dat");
		//else strcat(miscloop,"dnadynalignmiscloop.dat");

		if (opendat(loop, stackf, tstackh, tstacki,
			tloop, miscloop, danglef, int22,
			int21,coax, tstackcoax,
			coaxstack, tstack, tstackm, triloop,
			int11,hexaloop,tstacki23, tstacki1n, &data)) OK = true;
		else OK = false;

		pMainFrame = pMainFrame1;

		startpath = startpath1;

		singleinsert = Singleinsert;
		if (isRNA) {
			SetTitle("Dynalign");
		}
		else SetTitle ("Dynalign with DNA parameters");
	
		filenamedefined = false;

		T = 310.15; //Set the default temperature to 37 degrees C.
		Datapath = datapath;//Store the path to thermodynamic parameters in case the
							//user changes the temperature from 310.15 K.  (In which case
							//the enthalpy parameters need to be read.

}

void CDynDoc::namefile(char *Filename) {
	int i;

	filenamedefined = true;
	
	i = strlen(Filename);

	filename = new char [i+2];
	strcpy (filename,Filename);


}

int CDynDoc::newtemp() {
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

BOOL CDynDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;
	return TRUE;
}

CDynDoc::~CDynDoc()
{
	if(filenamedefined) delete[] filename;


}


BEGIN_MESSAGE_MAP(CDynDoc, CCTDoc)
	//{{AFX_MSG_MAP(CFoldDoc)
		// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CFoldDoc diagnostics

#ifdef _DEBUG
void CDynDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CDynDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CFoldDoc serialization

void CDynDoc::Serialize(CArchive& ar)
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
