// BiFoldDoc.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "BiFoldDoc.h"
#include "../src/defines.h"
#include "../src/algorithm.h"
#include "globals.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CBiFoldDoc

IMPLEMENT_DYNCREATE(CBiFoldDoc, CCTDoc)

CBiFoldDoc::CBiFoldDoc(bool isRNA, char *datapath, CWinApp *pMainFrame1, char *startpath1, bool fold, bool *Savefile)
{
	char loop[maxfil],stackf[maxfil], tstackh[maxfil], tstacki[maxfil],
		tloop[maxfil], miscloop[maxfil], danglef[maxfil], int22[maxfil],
		int21[maxfil],coax[maxfil], tstackcoax[maxfil],coaxstack[maxfil],
		tstack[maxfil], tstackm[maxfil],triloop[maxfil],int11[maxfil],hexaloop[maxfil],
		tstacki23[maxfil], tstacki1n[maxfil];

		ISRNA = isRNA;
	 
		getdat(loop, stackf, tstackh, tstacki,
			tloop, miscloop, danglef, int22,
			int21,coax, tstackcoax,
			coaxstack, tstack, tstackm, triloop,
			int11, hexaloop, tstacki23, tstacki1n, datapath, isRNA);

		if (opendat(loop, stackf, tstackh, tstacki,
			tloop, miscloop, danglef, int22,
			int21,coax, tstackcoax,
			coaxstack, tstack, tstackm, triloop,
			int11,hexaloop, tstacki23, tstacki1n,
			&data)) OK = true;
		else OK = false;

		pMainFrame = pMainFrame1;

		startpath = startpath1;


		if (isRNA&&fold) {
			SetTitle("RNA Fold");
		}
		else if (fold) SetTitle("DNA Fold");
		else if (isRNA) SetTitle("RNA Efn2");
		else SetTitle("DNA Efn2");

		savefile = Savefile;

		T = 310.15;//default folding temperature 
		Datapath = datapath; //keep track of the location of the thermodynamic parameters in case they are needed for loading enthalpy tables later
}

BOOL CBiFoldDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;
	return TRUE;
}

CBiFoldDoc::~CBiFoldDoc()
{
}

int CBiFoldDoc::newtemp() {
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


BEGIN_MESSAGE_MAP(CBiFoldDoc, CCTDoc)
	//{{AFX_MSG_MAP(CBiFoldDoc)
		// NOTE - the ClassWizard will add and remove mapping macros here.
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CBiFoldDoc diagnostics

#ifdef _DEBUG
void CBiFoldDoc::AssertValid() const
{
	CCTDoc::AssertValid();
}

void CBiFoldDoc::Dump(CDumpContext& dc) const
{
	CCTDoc::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CBiFoldDoc serialization

void CBiFoldDoc::Serialize(CArchive& ar)
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
// CBiFoldDoc commands
