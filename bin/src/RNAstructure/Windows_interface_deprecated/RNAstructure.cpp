// RNAstructure.cpp : Defines the class behaviors for the application.
//

//Copyright 2000,  David H. Mathews

#include "stdafx.h"
#include <WinDef.h>
#include <Windows.h>
#include <afxsettingsstore.h>
#include "RNAstructure.h"


#include "ChildFrm.h"
#include "RNAstructureDoc.h"
#include "RNAstructureView.h"
#include "FOLDVIEW.h"
#include "BiFoldView.h"
#include "CTDoc.h"
#include "Refolddoc.h"

#include "oligowalk.h"
#include "tprogressdialog.h"

#include "MyMDIChildWin.h"
#include "globals.h"
#include "drawview.h"
#include "drawdoc.h"
#include "drawmdichildwin.h"
#include "sequencedoc.h"
#include "oligoobject.h"
#include "oligoview.h"
#include "refoldview.h"
#include "dynrefoldview.h"
#include "dynaligndotplot.h"
#include "dotplot.h"
#include "dpdoc.h"
#include "dotplotview.h"

#include "pfdpdoc.h"
#include "boxplotview.h"
#include "pfformview.h"
#include "bimolpfview.h"
#include "allfoldview.h"
#include "RemovePseudo.h"

#include <direct.h>
#include <afxwin.h>
#include "SequenceEditor.h"
#include "Efn2_View.h"

#include "splash.h"
//#include "mixmatch.h"
#include <math.h>

#include "oligoscreen.h"
#include "oligoscreendoc.h"
#include "oligowalkrnaiview.h"


#include "Dynaligninterface.h"
#include "DynDoc.h"

#include "../src/outputconstraints.h"

#include "colorkey.h"
#include "colorkeydoc.h"
#include "SHAPEcolorkey.h"

#include "DynDotPlot.h"

#include "StochasticView.h"
#include "MEADialog.h"
#include "ProbKnotView.h"

#include "MultilignView.h"
#include "MultilignDoc.h"

#include "TurboFoldView.h"
#include "TurboFoldDoc.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif


#include <iostream>
using namespace std;

/////////////////////////////////////////////////////////////////////////////
// CRNAstructureApp

BEGIN_MESSAGE_MAP(CRNAstructureApp, CWinApp)
	
	//{{AFX_MSG_MAP(CRNAstructureApp)
	ON_COMMAND(ID_APP_ABOUT, OnAppAbout)
	ON_COMMAND(ID_FILE_FOLDDNASINGLESTRAND, OnFoldDNAsinglestrand)
	ON_COMMAND(ID_FILE_FOLDDNABIMOLECULAR, OnFoldDNAbimolecular)
	ON_COMMAND(ID_FILE_FOLDRNABIMOLECULAR, OnFoldRNAbimolecular)
	ON_COMMAND(ID_FILE_NEW, OnFileNew)
	ON_COMMAND(ID_FILE_EFN2, OnFileEfn2)
	ON_COMMAND(ID_FILE_EFN2DNA, OnFileEfn2dna)
	ON_COMMAND(ID_FILE_OPEN, OnFileOpen)
	ON_COMMAND(ID_FILE_OLIGOWALK, OnFileOligowalk)
	ON_COMMAND(ID_FILE_DYNALIGN, OnFileDynalign)
	//ON_COMMAND(ID_FILE_MIXMATCH, OnFileMixmatch)
	ON_COMMAND(ID_HELP_OPEN, OnHelpOpen)
	ON_COMMAND(ID_FILE_REFOLDFROMSAVEFILE, OnFileRefoldfromsavefile)
	ON_COMMAND(ID_FILE_DOTPLOTDYNALIGN,OnDynalignDotplot)
	ON_COMMAND(32802,OnDynalignRefold)
	//}}AFX_MSG_MAP
	// Standard file based document commands
	ON_COMMAND(ID_HELP_INDEX, OnHelpOpen)
	ON_COMMAND(ID_FILE_NEW, CWinApp::OnFileNew)
	ON_COMMAND(ID_FILE_OPEN, CWinApp::OnFileOpen)
	ON_COMMAND(ID_RNASINGLESTRAND, Rnasinglestrand)
	ON_COMMAND(ID_FILE_DRAW, OnDraw)
	// Standard print setup command
	ON_COMMAND(ID_FILE_PRINT_SETUP, CWinApp::OnFilePrintSetup)
	ON_COMMAND(ID_REFOLDFROMSINGLESTRAND, OnDynalignRefold)
	ON_COMMAND(ID_FILE_DOTPLOT, OnFileDotplot)
	ON_COMMAND(ID_FILE_DOTPLOTPARTITIONFUNCTION, OnFileDotplotpartitionfunction)
	ON_COMMAND(ID_FILE_PARTITIONFUNCTIONRNA, OnFilePartitionfunctionrna)
	
	ON_COMMAND(ID_FILE_GENERATEALLSUBOPTIMALSTRUCTURES, OnFileGenerateallsuboptimalstructures)
	ON_COMMAND(ID_FILE_OLIGOSCREEN, OnFileOligoscreen)
	ON_COMMAND(ID_FILE_OLIGOWALKFORSIRNA, OnFileOligowalkforsirna)
	ON_COMMAND(ID_FILE_PARTITIONFUNCTIONRNABIMOLECULAR, OnFilePartitionfunctionrnabimolecular)
	ON_COMMAND(ID_FILE_STOCHASTICTRACEBACK, OnFileStochastictraceback)
	ON_COMMAND(ID_FILE_PARTITIONFUNCTIONDNA, OnFilePartitionfunctiondna)
	ON_COMMAND(ID_FILE_DOTPLOTFROMTEXTFILE, OnFileDotplotfromtextfile)
	ON_COMMAND(ID_FILE_BREAKRNAPSEUDOKNOTS, OnFileBreakrnapseudoknots)
	ON_COMMAND(ID_RNA_PREDICTRNAMEASTRUCTURE, OnRNAMEA)
	ON_COMMAND(ID_DNA_PARTITIONFUNCTIONDNABIMOLECULAR, &CRNAstructureApp::OnDnaPartitionfunctiondnabimolecular)
	ON_COMMAND(ID_DNA_GENERATEALLDNASTRUCTURES, &CRNAstructureApp::OnDnaGeneratealldnastructures)
	ON_COMMAND(ID_DNA_STOCHASTICDNASAMPLING, &CRNAstructureApp::OnDnaStochasticdnasampling)
	ON_COMMAND(ID_DNA_DNADYNALIGN, &CRNAstructureApp::OnDnaDnadynalign)
	ON_COMMAND(ID_DNA_PREDICTDNAMEASTRUCTURE, &CRNAstructureApp::OnDnaPredictdnameastructure)
	ON_COMMAND(ID_RNA_PROBKNOT, &CRNAstructureApp::OnRnaProbknot)
	ON_COMMAND(ID_DNA_PROBKNOT, &CRNAstructureApp::OnDnaProbknot)
	ON_COMMAND(ID_RNA_RNAMULTILIGN, OnRNAMultilign)
	ON_COMMAND(ID_DNA_DNAMULTILIGN, &CRNAstructureApp::OnDnaDnamultilign)
	ON_COMMAND(ID_RNA_RNATURBOFOLD, &CRNAstructureApp::OnRnaRnaturbofold)
	ON_COMMAND(ID_DNA_DNATURBOFOLD, &CRNAstructureApp::OnDnaDnaturbofold)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CRNAstructureApp construction

CRNAstructureApp::CRNAstructureApp()
{
	// TODO: add construction code here,
	// Place all significant initialization in InitInstance
	
}

/////////////////////////////////////////////////////////////////////////////
// The one and only CRNAstructureApp object

CRNAstructureApp theApp;

/////////////////////////////////////////////////////////////////////////////
// CRNAstructureApp initialization

BOOL CRNAstructureApp::InitInstance()
{

	
	

	char dna[10];
	char conc[10];
	char *stopstring;

	{	// BLOCK: doc template registration
		// Register the document template.  Document templates serve
		// as the connection between documents, frame windows and views.
		// Attach this form to another document or frame window by changing
		// the document or frame class in the constructor below.
		CMultiDocTemplate* pNewDocTemplate = new CMultiDocTemplate(
			IDR_EFN2_VIEW_TMPL,
			RUNTIME_CLASS(CCTDoc),		// document class
			RUNTIME_CLASS(CMDIChildWnd),		// frame class
			RUNTIME_CLASS(CEfn2_View));		// view class
		AddDocTemplate(pNewDocTemplate);
	}


	{	// BLOCK: doc template registration
		// Register the document template.  Document templates serve
		// as the connection between documents, frame windows and views.
		// Attach this form to another document or frame window by changing
		// the document or frame class in the constructor below.
		CMultiDocTemplate* pNewDocTemplate = new CMultiDocTemplate(
			IDR_SEQUENCE,
			RUNTIME_CLASS(CCTDoc),		// document class
			RUNTIME_CLASS(CMDIChildWnd),		// frame class
			RUNTIME_CLASS(CSequenceEditor));		// view class
		AddDocTemplate(pNewDocTemplate);
	}

	char render[20];

	{	// BLOCK: doc template registration
		// Register the document template.  Document templates serve
		// as the connection between documents, frame windows and views.
		// Attach this form to another document or frame window by changing
		// the document or frame class in the constructor below.
		pFoldTemplate = new CMultiDocTemplate(
			IDR_FOLDVIEW_TMPL,
			RUNTIME_CLASS(CCTDoc),		// document class
			RUNTIME_CLASS(CMyMDIChildWin),		// frame class
			RUNTIME_CLASS(CFoldView));		// view class
		AddDocTemplate(pFoldTemplate);

		pAllFoldTemplate = new CMultiDocTemplate(
			IDR_ALLFOLDVIEW,
			RUNTIME_CLASS(CCTDoc),		// document class
			RUNTIME_CLASS(CMyMDIChildWin),		// frame class
			RUNTIME_CLASS(AllFoldView));		// view class
		AddDocTemplate(pAllFoldTemplate);

		pBreakPseudoTemplate = new CMultiDocTemplate(
			IDR_REMOVEPSEUDOFRAME,
			RUNTIME_CLASS(CCTDoc),		// document class
			RUNTIME_CLASS(CMyMDIChildWin),		// frame class
			RUNTIME_CLASS(CRemovePseudo));		// view class
		AddDocTemplate(pBreakPseudoTemplate);

		pReFoldTemplate = new CMultiDocTemplate(
			IDR_REFOLDVIEW_TMPL,
			RUNTIME_CLASS(CReFoldDoc),		// document class
			RUNTIME_CLASS(CMyMDIChildWin),		// frame class
			RUNTIME_CLASS(CRefoldView));		// view class
		AddDocTemplate(pReFoldTemplate);

		pBiFoldTemplate = new CMultiDocTemplate(
			IDR_BIFOLD,
			RUNTIME_CLASS(CBiFoldDoc),		// document class
			RUNTIME_CLASS(CMyMDIChildWin),		// frame class
			RUNTIME_CLASS(CBiFoldView));		// view class
		AddDocTemplate(pBiFoldTemplate);

		pDrawTemplate = new CMultiDocTemplate(
			IDR_DRAW,
			RUNTIME_CLASS(CDrawDoc),		// document class
			RUNTIME_CLASS(CDrawMDIChildWin),		// frame class
			RUNTIME_CLASS(CDrawView));		// view class
		AddDocTemplate(pDrawTemplate);

		pSequenceTemplate = new CMultiDocTemplate(
			IDR_SEQUENCE,
			RUNTIME_CLASS(CSequenceDoc),		// document class
			RUNTIME_CLASS(CMyMDIChildWin),		// frame class
			RUNTIME_CLASS(CSequenceEditor));		// view class
		AddDocTemplate(pSequenceTemplate);

		pEfn2Template = new CMultiDocTemplate(
			IDR_EFN2,
			RUNTIME_CLASS(CFoldDoc),
			RUNTIME_CLASS(CMyMDIChildWin),
			RUNTIME_CLASS(CEfn2_View));
		AddDocTemplate(pEfn2Template);

		pOligoTemplate = new CMultiDocTemplate(
			IDR_OLIGOWALK,
			RUNTIME_CLASS(COligoDoc),
			RUNTIME_CLASS(CMyMDIChildWin),
			RUNTIME_CLASS(COligoWalk));
		AddDocTemplate(pOligoTemplate);


		pOligoViewTemplate = new CMultiDocTemplate(
			IDR_OLIGOWALKVIEW,
			RUNTIME_CLASS(COligoDoc),
			RUNTIME_CLASS(CMDIChildWnd),
			RUNTIME_CLASS(COligoView));
		AddDocTemplate(pOligoViewTemplate);

		pOligosiRNATemplate = new CMultiDocTemplate(
			IDR_OLIGOWALK,
			RUNTIME_CLASS(COligoDoc),
			RUNTIME_CLASS(CMyMDIChildWin),
			RUNTIME_CLASS(COligoWalkRNAiView));
		AddDocTemplate(pOligosiRNATemplate);

		pOligoScreenTemplate = new CMultiDocTemplate(
			IDR_OLIGOSCREEN,
			RUNTIME_CLASS(COligoScreenDoc),
			RUNTIME_CLASS(CMyMDIChildWin),
			RUNTIME_CLASS(COligoScreen));
		AddDocTemplate(pOligoScreenTemplate);


		pDynalignTemplate = new CMultiDocTemplate(
			IDR_DYNALIGN,
			RUNTIME_CLASS(CDynDoc),
			RUNTIME_CLASS(CMyMDIChildWin),
			RUNTIME_CLASS(CDynalign));
			AddDocTemplate(pDynalignTemplate);

		//pMixmatchTemplate = new CMultiDocTemplate(
		//	IDR_MAINFRAME,
		//	RUNTIME_CLASS(CFoldDoc),
		//	RUNTIME_CLASS(CMyMDIChildWin),
		//	RUNTIME_CLASS(CMixMatch));
		//AddDocTemplate(pMixmatchTemplate);

		pDynAlignRefoldTemplate = new CMultiDocTemplate(
			IDR_DYNALIGN,
			RUNTIME_CLASS(CDynDoc),
			RUNTIME_CLASS(CMyMDIChildWin),
			RUNTIME_CLASS(DynRefoldView)
		);
		AddDocTemplate(pDynAlignRefoldTemplate);


		//partition function plot template:
		pPFDPTemplate = new CMultiDocTemplate(
			IDR_PLOT_PF,
			RUNTIME_CLASS(PFDPDoc),
			RUNTIME_CLASS(CDrawMDIChildWin),
			RUNTIME_CLASS(CBoxPlotView)
		);
		AddDocTemplate(pPFDPTemplate);

		pBimolPFTemplate = new CMultiDocTemplate(
			IDR_BIFOLDVIEW_TMPL,
			RUNTIME_CLASS(PFDPDoc),
			RUNTIME_CLASS(CDrawMDIChildWin),
			RUNTIME_CLASS(CBimolPFView)
		);
		AddDocTemplate(pBimolPFTemplate);

		

		pDPTemplate = new CMultiDocTemplate(
			IDR_PLOT,
			RUNTIME_CLASS(DPDoc),
			RUNTIME_CLASS(CDrawMDIChildWin),
			RUNTIME_CLASS(CDotPlotView)
		);
		AddDocTemplate(pDPTemplate);

		pPFTemplate = new CMultiDocTemplate(
			IDR_PFVIEW,
			RUNTIME_CLASS(CFoldDoc),
			RUNTIME_CLASS(CMyMDIChildWin),
			RUNTIME_CLASS(CPFFormView)
		);
		AddDocTemplate(pPFTemplate);
		
		pColorKeyTemplate= new CMultiDocTemplate(
			IDR_MAINFRAME,
			RUNTIME_CLASS(CColorKeyDoc),
			RUNTIME_CLASS(CMyMDIChildWin),
			RUNTIME_CLASS(CColorKey)
		);
		AddDocTemplate(pColorKeyTemplate);

		pSHAPEKeyTemplate= new CMultiDocTemplate(
			IDR_MAINFRAME,
			RUNTIME_CLASS(CColorKeyDoc),
			RUNTIME_CLASS(CMyMDIChildWin),
			RUNTIME_CLASS(CSHAPEKey)
		);
		AddDocTemplate(pSHAPEKeyTemplate);

		pStochasticTemplate= new CMultiDocTemplate(
			IDR_MAINFRAME,
			RUNTIME_CLASS(CFoldDoc),
			RUNTIME_CLASS(CMyMDIChildWin),
			RUNTIME_CLASS(CStochasticView)
		);
		AddDocTemplate(pStochasticTemplate);

		pMEATemplate= new CMultiDocTemplate(
			IDR_MAINFRAME,
			RUNTIME_CLASS(CFoldDoc),
			RUNTIME_CLASS(CMyMDIChildWin),
			RUNTIME_CLASS(CMEADialog)
		);
		AddDocTemplate(pMEATemplate);

		pPKTemplate= new CMultiDocTemplate(
			IDR_MAINFRAME,
			RUNTIME_CLASS(CFoldDoc),
			RUNTIME_CLASS(CMyMDIChildWin),
			RUNTIME_CLASS(CProbKnotView)
		);
		AddDocTemplate(pPKTemplate);

		pMultilignTemplate = new CMultiDocTemplate(
			IDR_EFN2,//For now, use efn2's menu because the only need is access to the temperature dialog
			RUNTIME_CLASS(CMultilignDoc),
			RUNTIME_CLASS(CMyMDIChildWin),
			RUNTIME_CLASS(MultilignView)
		);
		AddDocTemplate(pMultilignTemplate);


		pTurboFoldTemplate = new CMultiDocTemplate(
			IDR_EFN2,//For now, use efn2's menu because the only need is access to the temperature dialog
			RUNTIME_CLASS(CTurboFoldDoc),
			RUNTIME_CLASS(CMyMDIChildWin),
			RUNTIME_CLASS(CTurboFoldView)
		);
		AddDocTemplate(pTurboFoldTemplate);

	}


	CMultiDocTemplate* pDocTemplate;
	pDocTemplate = new CMultiDocTemplate(
		IDR_RNASTRTYPE,
		RUNTIME_CLASS(CCTDoc),
		RUNTIME_CLASS(CChildFrame), // custom MDI child frame
		RUNTIME_CLASS(CRNAstructureView));
	AddDocTemplate(pDocTemplate);


	AfxEnableControlContainer();

	// Standard initialization
	// If you are not using these features and wish to reduce the size
	//  of your final executable, you should remove from the following
	//  the specific initialization routines you do not need.

//#ifdef _AFXDLL
//	Enable3dControls();			// Call this when using MFC in a shared DLL
//#else
//	Enable3dControlsStatic();	// Call this when linking to MFC statically
//#endif

	// Change the registry key under which our settings are stored.
	// TODO: You should modify this string to be something appropriate
	// such as the name of your company or organization.
	//SetRegistryKey(_T("Local AppWizard-Generated Applications"));

	  // Load standard INI file options (including MRU)

	// Register the application's document templates.  Document templates
	//  serve as the connection between documents, frame windows and views.

	
	// create main MDI Frame window
	pMainFrame = new CMainFrame;
	if (!pMainFrame->LoadFrame(IDR_MAINFRAME))
		return FALSE;
	m_pMainWnd = pMainFrame;

	// Parse command line for standard shell commands, DDE, file open
	CCommandLineInfo cmdInfo;
	ParseCommandLine(cmdInfo);


	//commented out dhm 9/20/00 so that the program starts with no open documents
	// Dispatch commands specified on the command line
	//if (!ProcessShellCommand(cmdInfo))
	//	return FALSE;


	//custom code:
	//get the directory in which the program starts.  This is the location of the data files.
	


	//restore the path last used by the user
	/*SetRegistryKey(_T("RNAstructure"));
	strcpy(startpath,(GetProfileString("RNAstructure","startpath",datapath)).GetBuffer(10));
	strcpy(dna,(GetProfileString("RNAstructure","dna","Yes")).GetBuffer(1));

	if (!strcmp(dna,"Yes")) isdna = true;
	else isdna = false;

	strcpy(dna,(GetProfileString("RNAstructure","targetdna","Yes")).GetBuffer(1));

	if (!strcmp(dna,"Yes")) istargetdna = true;
	else istargetdna = false;
		
	option = (int) (GetProfileInt("RNAstructure","option",1)); 
	length = (int) (GetProfileInt("RNAstructure","length",18));
	strcpy(conc,(GetProfileString("RNAstructure","c","1e-6")).GetBuffer(1));


	c = strtod( conc, &stopstring);

	usesub = (int) (GetProfileInt("RNAstructure","usesub",1));



	

	strcpy(render,(GetProfileString("RNAstructure","render","clockwise")).GetBuffer(10));
	if (!strcmp("clockwise",render)) {

		clockwise = true;

	}
	else clockwise = false;

	//get a profile for dynalign
	strcpy(dna,(GetProfileString("RNAstructure","singleinsert","No")).GetBuffer(1));

	if (!strcmp(dna,"Yes")) singleinsert = true;
	else singleinsert = false;

	strcpy(dna,(GetProfileString("RNAstructure","savefile","No")).GetBuffer(1));

	if (!strcmp(dna,"Yes")) savefile = true;
	else savefile = false;

	
	strcpy(dna,(GetProfileString("RNAstructure","dynsavefile","No")).GetBuffer(1));

	if (!strcmp(dna,"Yes")) dynsavefile = true;
	else dynsavefile = false;*/



	//Restore Settings from Previous Uses of Program
	HKEY   hkey;
	DWORD  dwDisposition;
 
	DWORD dwType, dwSize;
	//int m_dwMaxFileSize = 16 * 1024;         // 16k
	
	// Set the default values in case nothing was previously set
	strcpy(startpath,datapath);
	isdna = false;
	istargetdna = false;
	option = 1; 
	length = 18;
	c = 1e-6;
	usesub = 1;
	clockwise = true;
	singleinsert = true;
	savefile = true;
	dynsavefile = true;

	_getcwd(datapath,_MAX_PATH);

	/*if(RegCreateKeyEx(HKEY_LOCAL_MACHINE, TEXT("Software\\University of Rochester\\RNAstructure"), 0, NULL, 0, 0, NULL, &hkey, &dwDisposition)== 
		ERROR_SUCCESS) {
	
		dwType = REG_SZ;
		dwSize = _MAX_PATH;
		RegQueryValueEx(hkey, TEXT("startpath"), NULL, &dwType, (PBYTE)&startpath, &dwSize);
        
		//dwType = REG_SZ;
		//dwSize = sizeof(m_szLastFileName);
		//RegQueryValueEx(hkey, TEXT("LastFileName"), NULL, &dwType, 
		//(PBYTE)&m_szLastFileName, &dwSize);
 
		RegCloseKey(hkey);
	}*/


	//Use the registry to read prior settings.
	//This is based on an example from: http://justins-fat-tire.blogspot.com/2009/10/mfc-reading-and-writing-from-registry.html
	
	// let's create a new pointer to our registry class 
	// Sadly the MSDN MFC documentation seems to be wrong
	// on the constructor logic.
	// This will put us under the HKEY_CURRENT_USER node
	// with full read/write rights.
	CSettingsStore* m_pRegistry = new CSettingsStore( FALSE, TRUE );
	CString string;
	int test;

	// First always check your pointer to make sure nothing
	// jinky has happened to it.
	// CreateKey opens the key or if it doesn't exist
	// creates a new key and sets it as the default current node.
	if( m_pRegistry!=NULL && m_pRegistry->Open( _T("Software\\Mathews_Lab\\RNAstructure") )!=0 ) {
		// read in the value from the node
		// again this will bring back a bool result for success or failure
		
		
		//m_pRegistry->Read( _T("test"), (int)test );

		m_pRegistry->Read(_T("startpath"), string);
		if (string!="") {
			//This means that startpath was recorded and all the other params should also have been recorded:

			strcpy(startpath,string.GetBuffer(0));
			m_pRegistry->Read(_T("isdna"), test);
			if (test>0) isdna = true;

			m_pRegistry->Read(_T("option"), option);
			
			m_pRegistry->Read(_T("istargetdna"), test);
			if (test>0) istargetdna = true;

			m_pRegistry->Read(_T("length"), length);
			
			m_pRegistry->Read(_T("usesub"), usesub);
			
			m_pRegistry->Read(_T("clockwise"), test);
			if (test>0) clockwise = true;

			m_pRegistry->Read(_T("singleinsert"), test);
			if (test>0) singleinsert = true;
			
			m_pRegistry->Read(_T("savefile"), test);
			if (test>0) savefile = true;

			m_pRegistry->Read(_T("dynsavefile"), test);
			if (test>0) dynsavefile = true;
			
			m_pRegistry->Read(_T("c"), string );
			c = strtod( string.GetBuffer(1), &stopstring);
			
	
	

		}

		//Added DATAPATH later, so be careful to make sure it is specified.
		m_pRegistry->Read(_T("DATAPATH"), string);
		if (string!="") {
			strcpy(datapath,string.GetBuffer(0));
			strcat(datapath,"\\");
		}

	
	else {
		
	}
		// always close our current connection when done
		m_pRegistry->Close();

	}
	



	if (__argc>1) {
		commandline();

	}
	else {
			// The main window has been initialized, so show and update it.
		pMainFrame->ShowWindow(m_nCmdShow);



		pMainFrame->ShowWindow(SW_SHOWMAXIMIZED);
		pMainFrame->UpdateWindow();


		//Show the splash dialog
		

		splash = new CSplash(pMainFrame);
		splash->Create(IDD_SPLASH, pMainFrame);
	}

	
	return TRUE;
}



int CRNAstructureApp::ExitInstance() 
{
		
	char conc[50];
	
	
	//write the startpath information to the registry:
	/*WriteProfileString("RNAstructure","startpath",startpath);//section,entry,value
	if (clockwise) WriteProfileString("RNAstructure","render","clockwise");
	else WriteProfileString("RNAstructure","render","counterclockwise");
	if (isdna) WriteProfileString("RNAstructure","dna","Yes");
	else WriteProfileString("RNAstructure","dna","No");
	if (istargetdna) WriteProfileString("RNAstructure","targetdna","Yes");
	else WriteProfileString("RNAstructure","targetdna","No");
	WriteProfileInt("RNAstructure","option",option);
	WriteProfileInt("RNAstructure","length",length);
	_gcvt(c,8,conc);
	WriteProfileString("RNAstructure","c",conc);
	WriteProfileInt("RNAstructure","usesub",usesub);
	if (singleinsert) WriteProfileString("RNAstructure","singleinsert","Yes");
	else WriteProfileString("RNAstructure","singleinsert","No");

	if (savefile) WriteProfileString("RNAstructure","savefile","Yes");
	else WriteProfileString("RNAstructure","savefile","No");

	if (dynsavefile) WriteProfileString("RNAstructure","dynsavefile","Yes");
	else WriteProfileString("RNAstructure","dynsavefile","No");*/


	//Write back to the registry:
	/*HKEY   hkey;
	DWORD  dwDisposition;
 
	DWORD dwType, dwSize;
 
	if(RegCreateKeyEx(HKEY_LOCAL_MACHINE, TEXT("Software\\University of Rochester\\RNAstructure"), 0, NULL, 0, 0, NULL, &hkey, &dwDisposition)== 
	ERROR_SUCCESS)
	{
		dwType = REG_SZ;
		dwSize = sizeof(startpath);
		RegSetValueEx(hkey, TEXT("startpath"), 0, dwType, (PBYTE)&startpath, dwSize);
 
		//dwType = REG_SZ;
		//dwSize = (_tcslen(m_szLastFileName) + 1) * sizeof(TCHAR);
		//RegSetValueEx(hkey, TEXT("LastFileName"), 0, dwType, 
		//(PBYTE)&m_szLastFileName, dwSize);
 
		RegCloseKey(hkey);
	}*/


	//Write to the registry; use example and comments from http://justins-fat-tire.blogspot.com/2009/10/mfc-reading-and-writing-from-registry.html

	// let's create a new pointer to our registry class 
	// Sadly the MSDN MFC documentation seems to be wrong
	// on the constructor logic.
	// This will put us under the HKEY_CURRENT_USER node
	// with full read/write rights.
	CSettingsStore* m_pRegistry = new CSettingsStore( FALSE, FALSE );

	
	CString string;

	// First always check your pointer to make sure nothing
	// jinky has happened to it.
	// CreateKey opens the key or if it doesn't exist
	// creates a new key and sets it as the default current node.
	if( m_pRegistry!=NULL && m_pRegistry->CreateKey( _T("Software\\Mathews_Lab\\RNAstructure") )!=0 ){
		
		// all the methods on the CSettingsClass return BOOL to
		// let you know if they succeeded
		// now we're under the node we want so we're going
		// to write some values, the Write methods looks for DWORD
		// so we need to do some casting to ensure our values get
		// written correctly
		
		//m_pRegistry->Write( _T("intValue"), (int)test );
		//m_pRegistry->Write( _T("boolValue"), (bool)testBool );

		string = startpath;
		m_pRegistry->Write( _T("startpath"), (LPCTSTR) string); 

		if(isdna) m_pRegistry->Write(_T("isdna"), (int) 1);
		else m_pRegistry->Write(_T("isdna"), (int) 0);

		m_pRegistry->Write(_T("option"), option);
			
		
		if(istargetdna) m_pRegistry->Write(_T("istargetdna"), (int) 1);
		else m_pRegistry->Write(_T("istargetdna"), (int) 0);

		m_pRegistry->Write(_T("length"), length);
			
		m_pRegistry->Write(_T("usesub"), usesub);
			
		if(clockwise) m_pRegistry->Write(_T("clockwise"), (int) 1);
		else m_pRegistry->Write(_T("clockwise"), (int) 0);

		if(singleinsert) m_pRegistry->Write(_T("singleinsert"), (int) 1);
		else m_pRegistry->Write(_T("singleinsert"), (int) 0);

		if(savefile) m_pRegistry->Write(_T("savefile"), (int) 1);
		else m_pRegistry->Write(_T("savefile"), (int) 0);

		if(dynsavefile) m_pRegistry->Write(_T("dynsavefile"), (int) 1);
		else m_pRegistry->Write(_T("dynsavefile"), (int) 0);

		_gcvt(c,8,conc);
		string = conc;
		m_pRegistry->Write(_T("c"), string );
		
	

		// always close our current connection when done
		m_pRegistry->Close();

	}



	//delete splash;

	return CWinApp::ExitInstance();
}


/////////////////////////////////////////////////////////////////////////////
// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	//{{AFX_DATA(CAboutDlg)
	enum { IDD = IDD_ABOUTBOX };
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CAboutDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	//{{AFX_MSG(CAboutDlg)
		// No message handlers
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:
	
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
	//{{AFX_DATA_INIT(CAboutDlg)
	//}}AFX_DATA_INIT
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CAboutDlg)
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
	//{{AFX_MSG_MAP(CAboutDlg)
		// No message handlers
	//}}AFX_MSG_MAP
	
END_MESSAGE_MAP()

// App command to run the dialog
void CRNAstructureApp::OnAppAbout()
{
	CAboutDlg aboutDlg;
	aboutDlg.DoModal();
}

/////////////////////////////////////////////////////////////////////////////
// CRNAstructureApp message handlers

//Fold RNA single strand
void CRNAstructureApp::Rnasinglestrand() {
	CFoldDoc *pFoldDocument = new CFoldDoc(true,datapath,this,startpath,true,&savefile);  //true indicates RNA parameters

	if (pFoldDocument->OK==true) openfoldwindow(pFoldDocument);
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;

}

void CRNAstructureApp::openfoldwindow(CFoldDoc *pFoldDocument) {
	/*CMultiDocTemplate* fold = new CMultiDocTemplate(
			IDR_FOLDVIEW_TMPL,
			RUNTIME_CLASS(CCTDoc),		// document class
			RUNTIME_CLASS(CMDIChildWnd),		// frame class
			RUNTIME_CLASS(CFoldview));		// view class*/
	
	//CMultiDocTemplate* fold;
	
	
	//CMDIChildWnd *pWin = new CMDIChildWnd;	

	CFrameWnd* pFrame = pFoldTemplate->CreateNewFrame(pFoldDocument,NULL);//NULL=pWin
	pFoldDocument->Frame=pFrame;
	pFoldDocument->menuframe = pMainFrame;
	

	pFoldTemplate->InitialUpdateFrame(pFrame,pFoldDocument);

	pFoldDocument->menuframe->GetMenu()->CheckMenuItem(ID_SUBFOLD,MF_UNCHECKED);

	//pWin->ModifyStyle(WS_THICKFRAME,0);//WS_MAXIMIZEBOX|,WS_OVERLAPPED
	


	return;
}

void CRNAstructureApp::openbifoldwindow(CBiFoldDoc *pBiFoldDocument) {
		

	CFrameWnd* pFrame = pBiFoldTemplate->CreateNewFrame(pBiFoldDocument,NULL);//NULL=pWin
	pBiFoldDocument->Frame=pFrame;
	pBiFoldDocument->menuframe = pMainFrame;
	pBiFoldDocument->menuframe->GetMenu()->CheckMenuItem(ID_FORCE_FORBIDUNIMOLECULARPAIRS,MF_UNCHECKED);

	pFoldTemplate->InitialUpdateFrame(pFrame,pBiFoldDocument);

	

	//pWin->ModifyStyle(WS_THICKFRAME,0);//WS_MAXIMIZEBOX|,WS_OVERLAPPED
	


	return;
}


void CRNAstructureApp::Draw(const char* ctfilename) {
	CDrawDoc *pDrawDocument;
	CRect rect,rect2;
	CPoint topleft,bottomright;
	int ht;



	pDrawDocument = new CDrawDoc(ctfilename,&clockwise,pMainFrame,this);
	CFrameWnd* pFrame = pDrawTemplate->CreateNewFrame(pDrawDocument,NULL);
	
	
	
	
	pMainFrame->GetClientRect( rect );
	pMainFrame->m_wndStatusBar.GetClientRect(rect2);
	
	//Here I do a pain-in-the-ass conversion to fix Microsoft's messed up code
	topleft = rect2.TopLeft();
	bottomright = rect2.BottomRight();
	ht = bottomright.y-topleft.y;
	rect.DeflateRect( 0,0,0,ht  );

	pMainFrame->m_wndToolBar.GetClientRect(rect2);

	topleft = rect2.TopLeft();
	bottomright = rect2.BottomRight();
	ht = bottomright.y-topleft.y+5;
	rect.DeflateRect( 0,0,5,ht  );


	
	pFrame->MoveWindow(rect);
	pDrawTemplate->InitialUpdateFrame(pFrame,pDrawDocument);

	CMenu* menu = pMainFrame->GetMenu( );
	
	
	if (*pDrawDocument->clockwise) menu->CheckMenuItem(ID_DRAW_CLOCKWISE,MF_CHECKED);
	else menu->CheckMenuItem(ID_DRAW_CLOCKWISE,MF_UNCHECKED);

	


}

void CRNAstructureApp::OnDraw() {
	//use a file open dialog box to get the name of a CT file to be drawn

	CFileDialog *filedialog;
	filedialog = new CFileDialog(TRUE,".ct",NULL,OFN_FILEMUSTEXIST,
		"CT Files (*.ct)|*.ct|All Files|*.*||");
	filedialog->m_ofn.lpstrInitialDir=startpath;
	if (filedialog->DoModal()==IDOK) {

		Draw((filedialog->GetPathName()).GetBuffer(0));
		//_getcwd(startpath,_MAX_PATH);

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
		strncpy(startpath,path.GetBuffer(1),i);
		*(startpath + i) ='\0';

	}
	delete filedialog;

	


}
void CRNAstructureApp::OnFoldDNAsinglestrand() 
{
	CFoldDoc *pFoldDocument = new CFoldDoc(false,datapath,this,startpath,true,&savefile);  //false indicates DNA parameters

	if (pFoldDocument->OK==true) openfoldwindow(pFoldDocument);
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;
	
}

void CRNAstructureApp::OnFoldDNAbimolecular() 
{
	CBiFoldDoc *pBiFoldDocument = new CBiFoldDoc(false,datapath,this,startpath,true,&savefile);  //false indicates DNA parameters

	if (pBiFoldDocument->OK==true) openbifoldwindow(pBiFoldDocument);
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Turner Lab website,\nhttp://rna.chem.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pBiFoldDocument;
	}
	return;
	
}

void CRNAstructureApp::OnFoldRNAbimolecular() 
{
	CBiFoldDoc *pBiFoldDocument = new CBiFoldDoc(true,datapath,this,startpath,true,&savefile);  //true indicates DNA parameters

	if (pBiFoldDocument->OK==true) openbifoldwindow(pBiFoldDocument);
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pBiFoldDocument;
	}
	return;
	
}

void CRNAstructureApp::OnFileNew() 
{
	
	CSequenceDoc *pSequenceDocument = new CSequenceDoc(pMainFrame,startpath,datapath,this);
	CFrameWnd* pFrame = pSequenceTemplate->CreateNewFrame(pSequenceDocument,NULL);//NULL=pWin
	

	pSequenceTemplate->InitialUpdateFrame(pFrame,pSequenceDocument);

	
}

void CRNAstructureApp::OnFileEfn2() 
{
	//Efn2 with RNA Rules
	CFoldDoc *pFoldDocument = new CFoldDoc(true,datapath,this,startpath,false); //true indicates RNA rules
	if (pFoldDocument->OK==true) openefn2window(pFoldDocument);
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;
}

void CRNAstructureApp::openefn2window(CFoldDoc *pFoldDocument) {
	CFrameWnd* pFrame = pEfn2Template->CreateNewFrame(pFoldDocument,NULL);//NULL=pWin
	pFoldDocument->Frame=pFrame;

	pEfn2Template->InitialUpdateFrame(pFrame,pFoldDocument);

}

void CRNAstructureApp::OnFileEfn2dna() 
{
	//Efn2 with DNA Rules
	CFoldDoc *pFoldDocument = new CFoldDoc(false,datapath,this,startpath,false); //false indicates DNA rules
	if (pFoldDocument->OK==true) openefn2window(pFoldDocument);
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;
	
}

void CRNAstructureApp::OnFileOpen() 
{
	// Open a sequence file
	CFileDialog *filedialog;
	filedialog = new CFileDialog(TRUE,".ct",NULL,OFN_FILEMUSTEXIST,
		"Sequence Files (*.seq)|*.seq|Plain Text Files (*.txt)|*.txt|Genbank Files (*.gen)|*.gen|All Files|*.*||");
	filedialog->m_ofn.lpstrInitialDir=startpath;

	if (filedialog->DoModal()==IDOK) {
		//_getcwd(startpath,_MAX_PATH);
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
		strncpy(startpath,path.GetBuffer(1),i);
		*(startpath + i) ='\0';
		CSequenceDoc *pSequenceDocument = new CSequenceDoc(pMainFrame,startpath,datapath,this,(filedialog->GetPathName()).GetBuffer(10));
		
		if (pSequenceDocument->check) {
			CFrameWnd* pFrame = pSequenceTemplate->CreateNewFrame(pSequenceDocument,NULL);//NULL=pWin
	

			pSequenceTemplate->InitialUpdateFrame(pFrame,pSequenceDocument);
		}
		else {
			delete pSequenceDocument;
			AfxMessageBox( "Sequence format was unrecognized!", 
				MB_OK|MB_ICONSTOP     );

		}
				

	}
	delete filedialog;

	

	
	
}

void CRNAstructureApp::FoldRNA(char *filename) {
	CFoldDoc *pFoldDocument = new CFoldDoc(true,datapath,this,startpath,true,&savefile);  //true indicates RNA parameters

	
	pFoldDocument->namefile(filename);
	
	
	if (pFoldDocument->OK==true) openfoldwindow(pFoldDocument);
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;



}

void CRNAstructureApp::FoldDNA(char *filename) {
	CFoldDoc *pFoldDocument = new CFoldDoc(false,datapath,this,startpath,true,&savefile);  //false indicates DNA parameters

	
	pFoldDocument->namefile(filename);
	
	
	if (pFoldDocument->OK==true) openfoldwindow(pFoldDocument);
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Turner Lab website,\nhttp://rna.chem.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;


}

void CRNAstructureApp::OnFileOligowalk() 
{

	COligoObject *oligoobject;

	oligoobject = new COligoObject(datapath,startpath,&isdna, &istargetdna, &option, 
		 &length,&c, &usesub);

	//Open the OligoWalk Dialog Box


		
	COligoDoc *pOligoDocument = new COligoDoc(oligoobject,pMainFrame,this); 
	CFrameWnd* pFrame = pOligoTemplate->CreateNewFrame(pOligoDocument,NULL);//NULL=pWin
	pOligoDocument->Frame=pFrame;
	pOligoDocument->pParent = this;
	pOligoDocument->SetTitle("OligoWalk");

	pOligoTemplate->InitialUpdateFrame(pFrame,pOligoDocument);
	


	


	
}


void CRNAstructureApp::OligoWalk(COligoObject *oligoobject) {
	
	COligoDoc *pOligoDocument;
	CString title;

	


	strcpy(startpath,oligoobject->startpath);
	isdna=oligoobject->isdna;
	istargetdna=oligoobject->istargetdna;
	option=oligoobject->option; 
	length=oligoobject->length;
	c=oligoobject->c;
	usesub=oligoobject->usesub;
	
	
	pOligoDocument = new COligoDoc(oligoobject,pMainFrame,this);
	CFrameWnd* pFrame = pOligoViewTemplate->CreateNewFrame(pOligoDocument,NULL);//NULL=pWin
	pOligoDocument->Frame=pFrame;
	title = "OligoWalk results for ";
	title+= oligoobject->ct.GetSequenceLabel().c_str();
	pOligoDocument->SetTitle(title);

	pOligoViewTemplate->InitialUpdateFrame(pFrame,pOligoDocument);

	CMenu* menu = pMainFrame->GetMenu( );
	
	
	menu->CheckMenuItem(ID_GRAPH_FREEENERGY_OVERALLANDDUPLEX,MF_CHECKED);
	

	
		


}


void CRNAstructureApp::DisplayOligoWalk(COligoObject *oligoobject) {
	


}


void CRNAstructureApp::OnFileDynalign() 
{
	 


	CDynDoc *pDynDocument = new CDynDoc(true,datapath,this,startpath,singleinsert);  //true indicates RNA parameters
	pDynDocument->checksave=&dynsavefile;

	if (pDynDocument->OK==true) {

		CFrameWnd* pFrame = pDynalignTemplate->CreateNewFrame(pDynDocument,NULL);//NULL=pWin
		pDynDocument->Frame=pFrame;

		pDynalignTemplate->InitialUpdateFrame(pFrame,pDynDocument);

	}
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pDynDocument;
	}
	return;



	
}

/*void CRNAstructureApp::OnFileMixmatch() 
{
	
	CFoldDoc *pFoldDocument = new CFoldDoc(true,datapath,this,startpath);  //true indicates RNA parameters

	if (pFoldDocument->OK==true) {
		CFrameWnd* pFrame = pMixmatchTemplate->CreateNewFrame(pFoldDocument,NULL);//NULL=pWin
		pFoldDocument->Frame=pFrame;

		pMixmatchTemplate->InitialUpdateFrame(pFrame,pFoldDocument);


	}
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;	


	
}*/

void CRNAstructureApp::commandline() {

	cout << "This graphical interface no longer supports command line parameters.  Please use the command line interfaces available from the Mathews lab.\n";

	/*structure ct,ct2,ct3;
	constrained = false;
	


	
	m_pMainWnd->ShowWindow(SW_HIDE);
	
	i = 2;

	if (!strcmp(__argv[1],"/fold")) {
		//fold a sequence
		isdna = false;
		seq1=false;
		bimol = false;
		
		number=20;
		percent = 20;
		window=-99;
		
		while (i<__argc) {
			if (!strcmp(__argv[i],"-s")) {
				if (!seq1) {
					strcpy(sequencefilename,__argv[i+1]);
					i = i+2;
					seq1= true;
				}
				else {
					strcpy(sequencefilename2,__argv[i+1]);
					i = i+2;
					bimol= true;


				}


			}
			else if (!strcmp(__argv[i],"-c")) {
				strcpy(ctfilename,__argv[i+1]);
				i = i+2;
			}
			
			else if (!strcmp(__argv[i],"-d")) {
				isdna=true;
				i++;

			}
			else if (!strcmp(__argv[i],"-n")) {
				
				number = atoi(__argv[i+1]);
				i=i+2;

			}
			else if (!strcmp(__argv[i],"-w")) {
				
				window = atoi(__argv[i+1]);
				i=i+2;

			}
			else if (!strcmp(__argv[i],"-p")) {
				
				percent = atoi(__argv[i+1]);
				i=i+2;

			}

			else if (!strcmp(__argv[i],"-con")) {

				strcpy(confilename,__argv[i+1]);
				i = i + 2;
				constrained = true;
			}
			

			


		}
		if (window==-99) {

				
			if (ct.numofbases>1200) {
   				window=20;
			}
			else if (ct.numofbases>800) {
   				window=15;
            
			}
			else if (ct.numofbases>500) {
   				window=11;
            
			}
			else if (ct.numofbases>300) {
   				window=7;
            
			}
			else if (ct.numofbases>120) {
   				window=5;
            
			}
			else if (ct.numofbases>50) {
   				window=3;
            
			}
			else window=2;



		}
		//Tell the user what is being done:
		cout << "\nRNAstructure structure prediction calculation with:\n";
		cout << "  Sequence: "<<sequencefilename<<"\n";
		if (bimol) cout << "  Sequence (2): "<<sequencefilename2<<"\n";
		cout << "  CT output to: " << ctfilename << "\n";
		cout << "  window: "<<window<<"\n";
		cout << "  maximum structures: "<<number<<"\n";
		cout << "  maximum % energy difference: "<<percent<<"\n";


		//now do the calculation:
		if (bimol) {
			openseq(&ct2,sequencefilename);
			openseq(&ct3,sequencefilename2);

			

			strcpy(ct.ctlabel[1],ct2.ctlabel[1]);
			//remove the new line at the end of ct.ctlabel[1]
			i = strlen(ct.ctlabel[1]);
			ct.ctlabel[1][i-1]='\0';

			strcat(ct.ctlabel[1],"_");
			strcat(ct.ctlabel[1],ct3.ctlabel[1]);

			
			ct.numofbases=ct2.numofbases+ct3.numofbases+3;
			ct.allocate(ct2.numofbases+ct3.numofbases+3);
			for (i=1;i<=ct2.numofbases;i++) {
				ct.numseq[i] = ct2.numseq[i];
				ct.nucs[i] = ct2.nucs[i];
				ct.hnumber[i] = ct2.hnumber[i];

			}
	
			for (i=1;i<=ct3.numofbases;i++) {
				ct.numseq[i+ct2.numofbases+3] = ct3.numseq[i];
				ct.nucs[i+ct2.numofbases+3] = ct3.nucs[i];
				ct.hnumber[i+ct2.numofbases+3] = ct3.hnumber[i];

			} 	
      
   
			ct.numseq[ct2.numofbases+1] = 5;
			ct.numseq[ct2.numofbases+2] = 5;
			ct.numseq[ct2.numofbases+3] = 5;

			ct.nucs[ct2.numofbases+1] = 'I';
			ct.nucs[ct2.numofbases+2] = 'I';
			ct.nucs[ct2.numofbases+3] = 'I';

			ct.hnumber[ct2.numofbases+1] = 0;
			ct.hnumber[ct2.numofbases+2] = 0;
			ct.hnumber[ct2.numofbases+3] = 0;


			ct.inter[0] = ct2.numofbases+1;
			ct.inter[1] = ct2.numofbases+2;
			ct.inter[2] = ct2.numofbases+3;

			ct.intermolecular = true;


		}
		else {
			openseq(&ct,sequencefilename);

		}
		
		
		
	
	 
		getdat(loop, stackf, tstackh, tstacki,
			tloop, miscloop, danglef, int22,
			int21,coax, tstackcoax,
			coaxstack, tstack, tstackm, triloop,
			int11, hexaloop,tstacki23, tstacki1n, datapath, !isdna);

		if (opendat(loop, stackf, tstackh, tstacki,
			tloop, miscloop, danglef, int22,
			int21,coax, tstackcoax,
			coaxstack, tstack, tstackm, triloop,
			int11,hexaloop,tstacki23, tstacki1n,&data)==0) {
		
			cout << "A Thermodynamic Data File was NOT Found.  Abort.\n";	
			return;

		}

		if (constrained) {
			//load the constraints into the ct file
			readconstraints(confilename,&ct);
		}


		dynamic(&ct,&data,number,percent,
			window);

		ctout (&ct,ctfilename);
		cout << "RNAstructure Calculation Complete\n";
		

	}
	else if (!strcmp(__argv[1],"/oligo")) {
		//oligo walk
		suboptimal = false;
		isdna=true;
		while (i<__argc) {
			if (!strcmp(__argv[i],"-c")) {
				strcpy(ctfilename,__argv[i+1]);
				i=i+2;

			}
			else if (!strcmp(__argv[i],"-r")) {
				strcpy(reportfilename,__argv[i+1]);
				i=i+2;

			}
			else if (!strcmp(__argv[i],"-m")) {
				mode = atoi(__argv[i+1]);
				i=i+2;

			}
			else if (!strcmp(__argv[i],"-l")) {
				length = atoi(__argv[i+1]);
				i=i+2;

			}
		
			else if (!strcmp(__argv[i],"/r")) {
				isdna=false;
				i++;

			}
			else if (!strcmp(__argv[i],"-s")) {
				suboptimal=true;
				i++;

			}
			else if (!strcmp(__argv[i],"-co")) {
				conc = atoi(__argv[i+1]);
				c = pow(1,((double) conc));
				i=i+2;

			}
		}

		cout << "\nRNAstructure OligoWalk calculation with:\n";
		cout << "  Target: "<<ctfilename<<"\n";
		cout << "  output to: " << reportfilename << "\n";
		if (isdna) cout << "  DNA Oligonucleotides\n";
		else cout << "  RNA Oligonucleotides\n";
		cout << "  mode: "<<mode<<"\n";
		
		openct(&ct,ctfilename);

		table = new int*[ct.numofbases - length + 2];

		for (i = 0; i < ct.numofbases - length + 2; i++) {
   			table[i] = new int[6];
		}

		if (isdna) {
			ddata = new datatable();
			hybriddata = new rddata;


		}
		thermo *helixstack = new thermo(datapath);
	


		getdat (loop2,stackf,tstackh,tstacki,tloop,miscloop,danglef,
   			int22,int21,coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,
			hexaloop,tstacki23, tstacki1n,datapath,true);

		if (opendat (loop2,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
   			coax,tstackcoax,coaxstack,tstack, tstackm, triloop, int11,hexaloop,
			tstacki23, tstacki1n, &data)==0) {

			cout << cout << "A Thermodynamic Data File was NOT Found.  Abort.\n";	
			return;


		}




		if (isdna) {
		

			getdat (loop2,stackf,tstackh,tstacki,tloop,miscloop,danglef,
   	      		int22,int21,coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,
				hexaloop,tstacki23, tstacki1n,datapath,false);

			if (opendat (loop2,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
   	      		coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,
				tstacki23, tstacki1n,ddata)==0) {


					cout << "A Thermodynamic Data File was NOT Found.  Abort.\n";	
					return;

		

			}

			
			strcpy(stackf,datapath);
			strcat (stackf,"\\");
			strcat (stackf,"stackdr.dat");
			readrd (hybriddata,stackf);

			strcpy(helixstack->DH,datapath);
			strcat(helixstack->DH,"\\");
			strcat(helixstack->DH,"stackdr.dh");

			strcpy(helixstack->DS,datapath);
			strcat(helixstack->DS,"\\");
			strcat(helixstack->DS,"stackdr.ds");

			strcpy(helixstack->HELIX,datapath);
			strcat(helixstack->HELIX,"\\");
			strcat(helixstack->HELIX,"helixdr.dat");
		
		}
	
		helixstack->read();

		olig(isdna, mode, &ct, length, 
			c, table, data, *ddata, 
			hybriddata,suboptimal,NULL,helixstack,
			1,ct.numofbases-length+1);


		report(reportfilename, &ct, table, 
			length,isdna, c, suboptimal);

		for (i = 0; i < ct.numofbases - length + 2; i++) {
   			delete[] table[i];
		}
		delete[] table;

		
		
		if (isdna) {
			delete ddata;
			//delete hybriddata;
			
			
		
		
		}
		delete helixstack;

		cout << "RNAstructure Calculation Complete\n";

	}
	else {
		//The first parameter was wrong, explain that to user
		cout <<"Usage RNAstructure [/type of calculation] [-modifiers]\nSee Online Help for more Information.";

	}
	
*/
	cout << "This version of the RNAstructure Windows User Interface no longer supports text-based commands.  Please use the new text-based interfaces available at http://rna.urmc.rochester.edu ./n";

	m_pMainWnd->SendMessage (WM_CLOSE);
}

void CRNAstructureApp::OnHelpOpen() 
{

	CString helpname;

	//helpname = datapath;
	//helpname+="//rnastructure.hlp";

//	::WinHelp(   pMainFrame->m_hWnd 
//,    helpname,    HELP_CONTENTS,
//    0);

	helpname = datapath;
	helpname+="\\html\\Contents.html";

	::ShellExecute(NULL,"open",helpname,NULL,NULL,SW_SHOWNORMAL);

	
	
}

void CRNAstructureApp::OnFileRefoldfromsavefile() 
{
	CReFoldDoc *pFoldDocument = new CReFoldDoc();
	//Open the refold window:
	CFrameWnd* pFrame = pReFoldTemplate->CreateNewFrame(pFoldDocument,NULL);//NULL=pWin
	pFoldDocument->Frame=pFrame;
	pFoldDocument->SetTitle("Refold from Save File");
	pFoldDocument->startingpath = startpath;
	pFoldDocument->pMainFrame = this;

	pFoldTemplate->InitialUpdateFrame(pFrame,pFoldDocument);



	
}

void CRNAstructureApp::OnDynalignRefold()

{
	CDynDoc *pDynDocument = new CDynDoc(true,datapath,this,startpath,singleinsert);  //true indicates RNA parameters
	

	

	CFrameWnd* pFrame = pDynAlignRefoldTemplate->CreateNewFrame(pDynDocument,NULL);//NULL=pWin
	pDynDocument->Frame=pFrame;

	pDynalignTemplate->InitialUpdateFrame(pFrame,pDynDocument);

	
	return;


}

// Create a Dynalign Dotplot
void CRNAstructureApp::OnDynalignDotplot(void)
{

	//The old code for the simple utility:
	/*
	//Open the DotPlot Dialog:
	CDynalignDotPlot *dotplot;

	dotplot = new CDynalignDotPlot();

	dotplot->Create(IDD_DYNALIGNDOTPLOT_DIALOG);*/

	//The new code for the GUI:

	//start by fetching the save file:
	CFileDialog *filedialog;
	filedialog = new CFileDialog(TRUE,".dsv",NULL,OFN_FILEMUSTEXIST,
		"Dynalign Save Files (*.dsv)|*.dsv|All Files|*.*||");
	filedialog->m_ofn.lpstrInitialDir=startpath;
	if (filedialog->DoModal()==IDOK) {

		DynDotPlot((filedialog->GetPathName()));
		//_getcwd(startpath,_MAX_PATH);	
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
		strncpy(startpath,path.GetBuffer(1),i);
		*(startpath + i) ='\0';

	}
	delete filedialog;



}

void CRNAstructureApp::DynDotPlot(CString Filename) {
	//Draw two dotplot windows, one for each sequence in a Dynalign calculation

		//open a dot plot view
	CString filename;
	CDynDotPlot *dpdoc1,*dpdoc2;
	CRect rect,rect2;
	CPoint topleft,bottomright;
	int ht;

	filename = Filename;

	dpdoc1 = new CDynDotPlot(filename,1,this);

	CFrameWnd* pFrame = pDPTemplate->CreateNewFrame(dpdoc1,NULL);

	pMainFrame->GetClientRect( rect );
	pMainFrame->m_wndStatusBar.GetClientRect(rect2);
	
	
	topleft = rect2.TopLeft();
	bottomright = rect2.BottomRight();
	ht = bottomright.y-topleft.y;
	rect.DeflateRect( 0,0,0,ht  );

	pMainFrame->m_wndToolBar.GetClientRect(rect2);

	topleft = rect2.TopLeft();
	bottomright = rect2.BottomRight();
	ht = bottomright.y-topleft.y+5;
	rect.DeflateRect( 0,0,5,ht  );

	pFrame->MoveWindow(rect);

	pDPTemplate->InitialUpdateFrame(pFrame,dpdoc1);

	////
	dpdoc2 = new CDynDotPlot(filename,2,this);

	CFrameWnd* pFrame2 = pDPTemplate->CreateNewFrame(dpdoc2,NULL);

	pFrame2->MoveWindow(rect);

	pDPTemplate->InitialUpdateFrame(pFrame2,dpdoc2);



}

void CRNAstructureApp::OnFileDotplot()
{
	//Open the DotPlot Dialog:
	/*CDotPlot *dotplot;//this was code for an old utility

	dotplot = new CDotPlot();

	dotplot->Create(IDD_DOTPLOT_DIALOG);*/

	CFileDialog *filedialog;
	filedialog = new CFileDialog(TRUE,".sav",NULL,OFN_FILEMUSTEXIST,
		"Save Files (*.sav)|*.sav|All Files|*.*||");
	filedialog->m_ofn.lpstrInitialDir=startpath;
	if (filedialog->DoModal()==IDOK) {

		DotPlot((filedialog->GetPathName()));
		//_getcwd(startpath,_MAX_PATH);	
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
		strncpy(startpath,path.GetBuffer(1),i);
		*(startpath + i) ='\0';

	}
	delete filedialog;
}

void CRNAstructureApp::OnFileDotplotpartitionfunction()
{

	//use a file open dialog box to get the name of a dot plot file to be drawn

	CFileDialog *filedialog;
	filedialog = new CFileDialog(TRUE,".pfs",NULL,OFN_FILEMUSTEXIST,
		"Partition Function Save Files (*.pfs)|*.pfs|All Files|*.*||");
	filedialog->m_ofn.lpstrInitialDir=startpath;
	if (filedialog->DoModal()==IDOK) {

		BoxPlot((filedialog->GetPathName()));
		//_getcwd(startpath,_MAX_PATH);	
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
		strncpy(startpath,path.GetBuffer(1),i);
		*(startpath + i) ='\0';

	}
	delete filedialog;

}

void CRNAstructureApp::BoxPlot(CString Filename) {
	//plot the dot plot


	//open a box plot view
	CString filename;
	PFDPDoc *pfdpdoc;
	CRect rect,rect2;
	CPoint topleft,bottomright;
	int ht;

	filename = Filename;

	pfdpdoc = new PFDPDoc(filename,this);

	//Check that the correct version of a file was read:
	if (pfdpdoc->arrayvalues == NULL) {

		//the file version doesn't match the current version
		AfxMessageBox( "Error: This save file was created with a different version of RNAstructure.  \nPlease redo the calculation.", 
		MB_OK|MB_ICONEXCLAMATION);

		delete pfdpdoc;

		return;

	}

	CFrameWnd* pFrame = pPFDPTemplate->CreateNewFrame(pfdpdoc,NULL);

	pMainFrame->GetClientRect( rect );
	pMainFrame->m_wndStatusBar.GetClientRect(rect2);
	
	
	topleft = rect2.TopLeft();
	bottomright = rect2.BottomRight();
	ht = bottomright.y-topleft.y;
	rect.DeflateRect( 0,0,0,ht  );

	pMainFrame->m_wndToolBar.GetClientRect(rect2);

	topleft = rect2.TopLeft();
	bottomright = rect2.BottomRight();
	ht = bottomright.y-topleft.y+5;
	rect.DeflateRect( 0,0,5,ht  );

	pFrame->MoveWindow(rect);

	pPFDPTemplate->InitialUpdateFrame(pFrame,pfdpdoc);




}

void CRNAstructureApp::DotPlot(CString Filename) {
	//plot the dot plot


	//open a dot plot view
	CString filename;
	DPDoc *dpdoc;
	CRect rect,rect2;
	CPoint topleft,bottomright;
	int ht;

	filename = Filename;

	dpdoc = new DPDoc(filename,this);

	CFrameWnd* pFrame = pDPTemplate->CreateNewFrame(dpdoc,NULL);

	pMainFrame->GetClientRect( rect );
	pMainFrame->m_wndStatusBar.GetClientRect(rect2);
	
	
	topleft = rect2.TopLeft();
	bottomright = rect2.BottomRight();
	ht = bottomright.y-topleft.y;
	rect.DeflateRect( 0,0,0,ht  );

	pMainFrame->m_wndToolBar.GetClientRect(rect2);

	topleft = rect2.TopLeft();
	bottomright = rect2.BottomRight();
	ht = bottomright.y-topleft.y+5;
	rect.DeflateRect( 0,0,5,ht  );

	pFrame->MoveWindow(rect);

	pDPTemplate->InitialUpdateFrame(pFrame,dpdoc);




}

void CRNAstructureApp::OnFilePartitionfunctionrna()
{
	//Do a partition function calculation for RNA

	CFoldDoc *pFoldDocument = new CFoldDoc(true,datapath,this,startpath,true,&savefile,true);  //true indicates RNA parameters

	
	
	if (pFoldDocument->OK==true) openpfwindow(pFoldDocument);
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;
}

void CRNAstructureApp::openpfwindow(CFoldDoc *pFoldDocument) {
	

	CFrameWnd* pFrame = pPFTemplate->CreateNewFrame(pFoldDocument,NULL);//NULL=pWin
	pFoldDocument->Frame=pFrame;

	pFoldDocument->menuframe = pMainFrame;
	

	pPFTemplate->InitialUpdateFrame(pFrame,pFoldDocument);

	//Uncheck the subfold menu item:
	pFoldDocument->menuframe->GetMenu()->CheckMenuItem(ID_SUBFOLD,MF_UNCHECKED);

}




void CRNAstructureApp::OnFileGenerateallsuboptimalstructures()
{


	CFoldDoc *pFoldDocument = new CFoldDoc(true,datapath,this,startpath,true,&savefile);  //true indicates RNA parameters

	if (pFoldDocument->OK==true) {
		CFrameWnd* pFrame = pAllFoldTemplate->CreateNewFrame(pFoldDocument,NULL);//NULL=pWin
		pFoldDocument->Frame=pFrame;

		pAllFoldTemplate->InitialUpdateFrame(pFrame,pFoldDocument);

	}
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;
	
}

void CRNAstructureApp::DisplayColorKeyWindow() {
	CColorKeyDoc *pCDocument = new CColorKeyDoc ();  

	
	CFrameWnd* pFrame = pColorKeyTemplate->CreateNewFrame(pCDocument,NULL);//NULL=pWin
	//pFoldDocument->Frame=pFrame;

	pColorKeyTemplate->InitialUpdateFrame(pFrame,pCDocument);


}

void CRNAstructureApp::DisplaySHAPEKeyWindow() {
	CColorKeyDoc *pCDocument = new CColorKeyDoc ();  

	
	CFrameWnd* pFrame = pSHAPEKeyTemplate->CreateNewFrame(pCDocument,NULL);//NULL=pWin
	//pFoldDocument->Frame=pFrame;

	pSHAPEKeyTemplate->InitialUpdateFrame(pFrame,pCDocument);


}

void CRNAstructureApp::OnFileOligoscreen()
{
	//Open the OligoScreen Window:
	COligoScreenDoc *pOSD = new COligoScreenDoc(datapath,startpath);  //true indicates RNA parameters
	

	CFrameWnd* pFrame = pOligoScreenTemplate->CreateNewFrame(pOSD,NULL);//NULL=pWin
	pOSD->Frame=pFrame;

	pOligoScreenTemplate->InitialUpdateFrame(pFrame,pOSD);


}

void CRNAstructureApp::OnFileOligowalkforsirna()
{
	//Open a version of OligoWalk geared to siRNA design
	COligoObject *oligoobject;

	isdna = false;//set isdna false because this is siRNA design -- using RNA
	oligoobject = new COligoObject(datapath,startpath,&isdna, &istargetdna, &option, 
		 &length,&c, &usesub);

	//Open the OligoWalk Dialog Box	
	COligoDoc *pOligoDocument = new COligoDoc(oligoobject,pMainFrame,this); 
	CFrameWnd* pFrame = pOligosiRNATemplate->CreateNewFrame(pOligoDocument,NULL);//NULL=pWin
	pOligoDocument->Frame=pFrame;
	pOligoDocument->pParent = this;
	pOligoDocument->SetTitle("OligoWalk for siRNA");
	oligoobject->siRNA = true; //set true to indicate this is an siRNA design

	pOligosiRNATemplate->InitialUpdateFrame(pFrame,pOligoDocument);
	

}

void CRNAstructureApp::OnFilePartitionfunctionrnabimolecular()
{
	// Bimolecular Partition Function:

	CFoldDoc *pFoldDocument = new CFoldDoc(true,datapath,this,startpath,true,&savefile,true);  //true indicates RNA parameters


	
	if (pFoldDocument->OK==true) openbimolpfwindow(pFoldDocument);
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;
}

void CRNAstructureApp::openbimolpfwindow(CFoldDoc *pFoldDocument) {
	

	CFrameWnd* pFrame = pBimolPFTemplate->CreateNewFrame(pFoldDocument,NULL);//NULL=pWin
	pFoldDocument->Frame=pFrame;

	

	pBimolPFTemplate->InitialUpdateFrame(pFrame,pFoldDocument);

}


void CRNAstructureApp::OnFileStochastictraceback()
{
	//Perform a stochastic traceback:

	CFoldDoc *pFoldDocument = new CFoldDoc(true,datapath,this,startpath,true,&savefile);  //false indicates DNA parameters

	if (pFoldDocument->OK==true) {
		CFrameWnd* pFrame = pStochasticTemplate->CreateNewFrame(pFoldDocument,NULL);//NULL=pWin
		pFoldDocument->Frame=pFrame;
		pFoldDocument->menuframe = pMainFrame;
		pStochasticTemplate->InitialUpdateFrame(pFrame,pFoldDocument);

	

	}
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;
	
	

}

void CRNAstructureApp::OnFilePartitionfunctiondna()
{
	//Partition Function for a DNA Sequence:
	CFoldDoc *pFoldDocument = new CFoldDoc(false,datapath,this,startpath,true,&savefile,true);  //false indicates DNA parameters


	

	if (pFoldDocument->OK==true) openpfwindow(pFoldDocument);
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;

}

void CRNAstructureApp::OnFileDotplotfromtextfile()
{
	// Build a dot plot based on a text file in the format i j energy


	CFileDialog *filedialog;
	filedialog = new CFileDialog(TRUE,".txt",NULL,OFN_FILEMUSTEXIST,
		"Dot Plot File (*.dp)|*.dp|Plain Text File (*.txt)|*.txt|All Files|*.*||");
	filedialog->m_ofn.lpstrInitialDir=startpath;
	if (filedialog->DoModal()==IDOK) {
	
		//open a dot plot view
		CString filename;
		DPDoc *dpdoc;
		CRect rect,rect2;
		CPoint topleft,bottomright;
		int ht;

		filename = filedialog->GetPathName();

		dpdoc = new DPDoc(filename,this,true);

		CFrameWnd* pFrame = pDPTemplate->CreateNewFrame(dpdoc,NULL);

		pMainFrame->GetClientRect( rect );
		pMainFrame->m_wndStatusBar.GetClientRect(rect2);
	
	
		topleft = rect2.TopLeft();
		bottomright = rect2.BottomRight();
		ht = bottomright.y-topleft.y;
		rect.DeflateRect( 0,0,0,ht  );

		pMainFrame->m_wndToolBar.GetClientRect(rect2);

		topleft = rect2.TopLeft();
		bottomright = rect2.BottomRight();
		ht = bottomright.y-topleft.y+5;
		rect.DeflateRect( 0,0,5,ht  );

		pFrame->MoveWindow(rect);

		pDPTemplate->InitialUpdateFrame(pFrame,dpdoc);
	}
	delete filedialog;
}

void CRNAstructureApp::OnFileBreakrnapseudoknots()
{
	CFoldDoc *pFoldDocument = new CFoldDoc(true,datapath,this,startpath,true,&savefile);  //true indicates RNA parameters

	if (pFoldDocument->OK==true) {
		CFrameWnd* pFrame = pBreakPseudoTemplate->CreateNewFrame(pFoldDocument,NULL);//NULL=pWin
		pFoldDocument->Frame=pFrame;

		pBreakPseudoTemplate->InitialUpdateFrame(pFrame,pFoldDocument);

	}
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;
	

}


void CRNAstructureApp::OnRNAMEA()
{
	//Perform a stochastic traceback:

	CFoldDoc *pFoldDocument = new CFoldDoc(true,datapath,this,startpath,true,&savefile);  //true indicates RNA parameters

	if (pFoldDocument->OK==true) {
		CFrameWnd* pFrame = pMEATemplate->CreateNewFrame(pFoldDocument,NULL);//NULL=pWin
		pFoldDocument->Frame=pFrame;
		pFoldDocument->menuframe = pMainFrame;
		pMEATemplate->InitialUpdateFrame(pFrame,pFoldDocument);

	

	}
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;
	
	

}

//Handle a request for DNA bimolecular partition functions
void CRNAstructureApp::OnDnaPartitionfunctiondnabimolecular()
{

	// Bimolecular Partition Function:

	CFoldDoc *pFoldDocument = new CFoldDoc(false,datapath,this,startpath,true,&savefile,true);  //first false indicates DNA

	
	if (pFoldDocument->OK==true) openbimolpfwindow(pFoldDocument);
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;
	
}

void CRNAstructureApp::OnDnaGeneratealldnastructures()
{
	//Predict All Suboptimal Structures Using DNA Parameters
	CFoldDoc *pFoldDocument = new CFoldDoc(false,datapath,this,startpath,true,&savefile);  //false indicates DNA parameters

	if (pFoldDocument->OK==true) {
		CFrameWnd* pFrame = pAllFoldTemplate->CreateNewFrame(pFoldDocument,NULL);//NULL=pWin
		pFoldDocument->Frame=pFrame;

		pAllFoldTemplate->InitialUpdateFrame(pFrame,pFoldDocument);

	}
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;
}

void CRNAstructureApp::OnDnaStochasticdnasampling()
{
	//Perform DNA Stochastic Sampling.  Note that the thermodynamic parameters are stored with the pfs file, so this requires the exact same calla s the RNA version.


	CFoldDoc *pFoldDocument = new CFoldDoc(false,datapath,this,startpath,true,&savefile);  //false indicates DNA parameters

	if (pFoldDocument->OK==true) {
		CFrameWnd* pFrame = pStochasticTemplate->CreateNewFrame(pFoldDocument,NULL);//NULL=pWin
		pFoldDocument->Frame=pFrame;
		pFoldDocument->menuframe = pMainFrame;
		pStochasticTemplate->InitialUpdateFrame(pFrame,pFoldDocument);

	

	}
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;
	
}

void CRNAstructureApp::OnDnaDnadynalign()
{
	// Perform Dynalign with DNA parameters:
	CDynDoc *pDynDocument = new CDynDoc(false,datapath,this,startpath,singleinsert);  //false indicates DNA parameters
	pDynDocument->checksave=&dynsavefile;

	if (pDynDocument->OK==true) {

		CFrameWnd* pFrame = pDynalignTemplate->CreateNewFrame(pDynDocument,NULL);//NULL=pWin
		pDynDocument->Frame=pFrame;

		pDynalignTemplate->InitialUpdateFrame(pFrame,pDynDocument);

	}
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pDynDocument;
	}
	return;
	

}

void CRNAstructureApp::OnDnaPredictdnameastructure()
{
	//These can be handled exactly the same as RNA because the save file has the thermodynamic parameters stored.
	

	CFoldDoc *pFoldDocument = new CFoldDoc(false,datapath,this,startpath,true,&savefile);  //false indicates DNA parameters

	if (pFoldDocument->OK==true) {
		CFrameWnd* pFrame = pMEATemplate->CreateNewFrame(pFoldDocument,NULL);//NULL=pWin
		pFoldDocument->Frame=pFrame;
		pFoldDocument->menuframe = pMainFrame;
		pMEATemplate->InitialUpdateFrame(pFrame,pFoldDocument);

	

	}
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;
}


void CRNAstructureApp::OnRnaProbknot() {

	//Use an RNA pFoldDocument just to get the correct name in the frame, otherwise the parameters are superfluous.
	CFoldDoc *pFoldDocument = new CFoldDoc(true,datapath,this,startpath,true,&savefile);  //false indicates DNA parameters

	if (pFoldDocument->OK==true) {
		CFrameWnd* pFrame = pPKTemplate->CreateNewFrame(pFoldDocument,NULL);//NULL=pWin
		pFoldDocument->Frame=pFrame;
		pFoldDocument->menuframe = pMainFrame;
		pPKTemplate->InitialUpdateFrame(pFrame,pFoldDocument);

	

	}
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;

}

//Perform a ProbKnot Calculation
void CRNAstructureApp::OnDnaProbknot() {

	//Use a DNA pFoldDocument just to get the correct name in the frame, otherwise the parameters are superfluous.
	CFoldDoc *pFoldDocument = new CFoldDoc(false,datapath,this,startpath,true,&savefile);  //false indicates DNA parameters

	if (pFoldDocument->OK==true) {
		CFrameWnd* pFrame = pPKTemplate->CreateNewFrame(pFoldDocument,NULL);//NULL=pWin
		pFoldDocument->Frame=pFrame;
		pFoldDocument->menuframe = pMainFrame;
		pPKTemplate->InitialUpdateFrame(pFrame,pFoldDocument);

	

	}
	else {
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
		delete pFoldDocument;
	}
	return;

}


//Perform Multilign with RNA Rules
void CRNAstructureApp::OnRNAMultilign() {


	
	CMultilignDoc *pmultidoc = new CMultilignDoc(true,datapath,this,startpath);  //true indicates RNA parameters
	

	

	CFrameWnd* pFrame = pMultilignTemplate->CreateNewFrame(pmultidoc,NULL);//NULL=pWin
	pmultidoc->Frame=pFrame;

	pMultilignTemplate->InitialUpdateFrame(pFrame,pmultidoc);

	
	return;
}

//Perform Multilign with DNA Rules
void CRNAstructureApp::OnDnaDnamultilign()
{
	CMultilignDoc *pmultidoc = new CMultilignDoc(false,datapath,this,startpath);  //false indicates DNA parameters
	

	

	CFrameWnd* pFrame = pMultilignTemplate->CreateNewFrame(pmultidoc,NULL);//NULL=pWin
	pmultidoc->Frame=pFrame;

	pMultilignTemplate->InitialUpdateFrame(pFrame,pmultidoc);

	
	return;
}



//Perform TurboFold using RNA rules:
void CRNAstructureApp::OnRnaRnaturbofold()
{
	CTurboFoldDoc *pturbodoc = new CTurboFoldDoc(true,datapath,this,startpath);  //true indicates RNA parameters
	

	

	CFrameWnd* pFrame = pTurboFoldTemplate->CreateNewFrame(pturbodoc,NULL);//NULL=pWin
	pturbodoc->Frame=pFrame;

	pTurboFoldTemplate->InitialUpdateFrame(pFrame,pturbodoc);

	
	return;
}


//Perform TurboFold using DNA rules:
void CRNAstructureApp::OnDnaDnaturbofold()
{
	CTurboFoldDoc *pturbodoc = new CTurboFoldDoc(false,datapath,this,startpath);  //false indicates DNA parameters
	

	

	CFrameWnd* pFrame = pTurboFoldTemplate->CreateNewFrame(pturbodoc,NULL);//NULL=pWin
	pturbodoc->Frame=pFrame;

	pTurboFoldTemplate->InitialUpdateFrame(pFrame,pturbodoc);

	
}
