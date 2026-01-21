// OligoScreen.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "OligoScreen.h"
#include "globals.h"
#include <direct.h>

#include "temp_Dialog.h"
#include "../src/OligoScreenCalc.h"



using namespace std;





//run oligoscreen in a seperate thread
UINT oligoscreen(LPVOID pParam) {//(char *infilename, char *outfilename, bool isRNA, char *path, TProgressDialog *progress) {
	
	
	int i,j,k,l,energy;
	datatable data;
	datatable *enthalpy;
	char loop[maxfil], stackf[maxfil], tstackh[maxfil], tstacki[maxfil],
		tloop[maxfil], miscloop[maxfil], danglef[maxfil], int22[maxfil],
		int21[maxfil], coax[maxfil], tstackcoax[maxfil],
		coaxstack[maxfil], tstack[maxfil], tstackm[maxfil], triloop[maxfil],
		int11[maxfil],hexaloop[maxfil],tstack23[maxfil],tstacki1n[maxfil];
	rddata *hybriddata,*enthalpyhybrid;

	COligoScreenObject* pObject = (COligoScreenObject*) pParam;
	COligoScreen* pView = (COligoScreen*) pObject->parent;

	//open the data files for folding
	getdat(loop,stackf, tstackh, tstacki,
		tloop, miscloop, danglef, int22,
		int21,coax, tstackcoax,
		coaxstack, tstack, tstackm, triloop,
		int11, hexaloop, tstack23, tstacki1n, pObject->path, pObject->isRNA);

	

	opendat(loop, stackf, tstackh, tstacki,
		tloop, miscloop, danglef, int22,
		int21,coax, tstackcoax,
		coaxstack, tstack, tstackm, triloop,
		int11,hexaloop,tstack23,tstacki1n,&data);

	//check to see if the temperature has been changed.
	if (pObject->T<310||pObject->T>311) {
		getenthalpydat(loop, stackf, tstackh, tstacki,
			tloop, miscloop, danglef, int22,
			int21,coax, tstackcoax,
			coaxstack, tstack, tstackm, triloop,
			int11, hexaloop, tstack23, tstacki1n, pObject->path, pObject->isRNA);

		enthalpy = new datatable();

		opendat(loop, stackf, tstackh, tstacki,
			tloop, miscloop, danglef, int22,
			int21,coax, tstackcoax,
			coaxstack, tstack, tstackm, triloop,
			int11,hexaloop,tstack23, tstacki1n, enthalpy);
	
		
		dG_T(pObject->T,data,*enthalpy,data);

		delete enthalpy;	
		

	}

	if (!pObject->isRNA) {
		//This is DNA oligos

		hybriddata = new rddata;
		strcpy(stackf,pObject->path);
		strcat(stackf,"\\");
		strcat (stackf,"stackdr.dat");
		readrd(hybriddata,stackf);
		
		if (pObject->T<310||pObject->T>311) {
			//The temperature has changed.
			//Read the enthalpy data into a rddata.
			
			strcpy(stackf,pObject->path);
			strcat(stackf,"\\");
			strcat (stackf,"stackdr.dh");
			enthalpyhybrid = new rddata;
			readrd(enthalpyhybrid,stackf);
      		

			for (i=0;i<5;i++) {
				for (j=0;j<5;j++) {
					for (k=0;k<5;k++) {
						for (l=0;l<5;l++) {
							hybriddata->stack[i][j][k][l]=Tscale(pObject->T,hybriddata->stack[i][j][k][l],enthalpyhybrid->stack[i][j][k][l]);
						}
					}
				}
			}
			delete enthalpyhybrid;
		}

	}
	else {
		//This is RNA oligos
		hybriddata= NULL;
	}

	//call the backend function
	OligoScreenCalc(pObject->infilename, pObject->outfilename, &data, hybriddata);


	if (!pObject->isRNA) 
		delete hybriddata;

	::PostMessage(pView->m_hWnd,ID_FOLDDONE,0,0);


	return 0;//thread completed successfully



}

// COligoScreen

IMPLEMENT_DYNCREATE(COligoScreen, CFormView)

COligoScreen::COligoScreen()
	: CFormView(COligoScreen::IDD)
	, m_listname(_T(""))
	, m_outname(_T(""))
{

	T = (float) 310.15;//The default temperature of the calculation.

}

COligoScreen::~COligoScreen()
{
}

void COligoScreen::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_OLIGOMERLISTNAME, m_listname);
	DDX_Text(pDX, IDC_OUTFILE_NAME, m_outname);
}

BEGIN_MESSAGE_MAP(COligoScreen, CFormView)
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
	ON_BN_CLICKED(IDC_LIST, OnBnClickedList)
	ON_BN_CLICKED(IDC_OUT, OnBnClickedOut)
	ON_MESSAGE(ID_FOLDDONE, DoneFolding)
	ON_COMMAND(ID_TEMPERATURE, OnTemperature)
END_MESSAGE_MAP()


// COligoScreen diagnostics

#ifdef _DEBUG
void COligoScreen::AssertValid() const
{
	CFormView::AssertValid();
}

void COligoScreen::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}
#endif //_DEBUG


// COligoScreen message handlers

void COligoScreen::OnInitialUpdate(void)
{
	CFormView::OnInitialUpdate();
	CButton *dna;

	
	CMenu* pSysMenu = GetSystemMenu(FALSE);
	


	dna = (CButton*) GetDlgItem( IDC_DNA );
	dna->SetCheck(1);


	ResizeParentToFit();
	
	

}

void COligoScreen::OnBnClickedOk()
{
		//Extract all the information and run OligoScreen
	bool isRNA;
	CButton *rna;
	


	rna = (CButton*) GetDlgItem( IDC_RNA );
	


	if (rna->GetCheck()) {
		isRNA=true;

	}
	else isRNA=false;

	if (m_listname!="") { 
		//make sure the user has chosen a list file 


		CRect *rect;
		rect = new CRect(10,40,210,120);
		progress = new TProgressDialog();
		progress->Create(NULL,"Folding the Oligos...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);


		oligoscreenobject.infilename = m_listname.GetBuffer(1);
		oligoscreenobject.outfilename = m_outname.GetBuffer(1);
		oligoscreenobject.isRNA = isRNA;
		oligoscreenobject.path=GetOligoScreenDocument()->datapath;
		oligoscreenobject.progress = progress;
		oligoscreenobject.parent=this;
		oligoscreenobject.T=T;

		//start a new thread for the background work:
		AfxBeginThread(oligoscreen,&oligoscreenobject);
			
		//if(oligoscreen(m_listname.GetBuffer(1), m_outname.GetBuffer(1),isRNA,GetOligoScreenDocument()->datapath,NULL)==0) 
		//	AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
		//	MB_OK|MB_ICONSTOP);
		
	}
	else return;
	
	
}

COligoScreenDoc *COligoScreen::GetOligoScreenDocument() {
	
	return ((COligoScreenDoc*) GetDocument());


}


void COligoScreen::OnBnClickedList()
{
		//Get the name of the list file:
	char *outname;
	short int i;
	CFileDialog *filedialog;
	
	

	
	
	filedialog = new CFileDialog(TRUE,0,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"List Files (*.lis)|*.lis|All Files|*.*||");

	
	
	if (filedialog->DoModal()==IDOK) {
		
		m_listname=(filedialog->GetPathName()).GetBuffer(30);
		i = m_listname.GetLength();
		
		outname = new char[i+4];//allocate enough space so that 
														//three characters can be added 
														//to the name if necessary
		strcpy(outname,m_listname.GetBuffer(10));
		//count the characters to the .
		
		while(i>=0){
			
			if (outname[i]=='.') break;
			i--;
		}
		if (i==0) i = m_listname.GetLength();
		strcpy(outname+i+1,"out\0");
		m_outname=outname;
		
		delete[] outname;
		UpdateData(FALSE);
		
		

	}
	delete filedialog;

}

void COligoScreen::OnBnClickedOut()
{
		//The user is specifying the CT file name explicitly
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".out",NULL,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"out Files (*.out)|*.out|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		m_outname=(filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);		

	}
	delete filedialog;
}

LRESULT COligoScreen::DoneFolding(WPARAM wParam, LPARAM lParam) {
	
	progress->SendMessage (WM_CLOSE);
	
	GetOligoScreenDocument()->Frame->SendMessage(WM_CLOSE);
	
	return 0;


}

void COligoScreen::OnTemperature()
{
	//Allow the user to specify a new temperature.
	CTemp_Dialog *temp;
	temp=new CTemp_Dialog(&T);

	temp->DoModal();
	delete temp;
	


}
