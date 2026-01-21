

// RNAstructure.h : main header file for the RNASTRUCTURE application
//

#if !defined(AFX_RNASTRUCTURE_H__56CD85E4_8D8A_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_RNASTRUCTURE_H__56CD85E4_8D8A_11D4_9F32_00C0F02A5F5D__INCLUDED_


#include "FoldDoc.h"
#include "MainFrm.h"
#include "BiFoldDoc.h"
#include "TProgressDialog.h"
#include "../src/structure.h"
#include "../src/rna_library.h"
#include "../src/intermolecular.h"
#include "oligoobject.h"



#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"       // main symbols
#include "splash.h"

/////////////////////////////////////////////////////////////////////////////
// CRNAstructureApp:
// See RNAstructure.cpp for the implementation of this class
//


//UINT OligoWalkProc( LPVOID pParam );






class CRNAstructureApp : public CWinApp
{
public:
	CRNAstructureApp();
	void Rnasinglestrand();
	void openfoldwindow(CFoldDoc *pFoldDocument);
	void openbifoldwindow(CBiFoldDoc *pBiFoldDocument); 
	void openefn2window(CFoldDoc *pFoldDocument);
	void DisplayOligoWalk(COligoObject *oligoobject);
	void OligoWalk(COligoObject *oligoobject); 
	char datapath[_MAX_PATH];//keep track of the startup path so that the program knows
							//where to fetch the data files
	void Draw(const char* filename);
	void OnDraw();
	char startpath[_MAX_PATH];//keep track of the path that the user last used to open a file
							//this is saved and restored from the registry
	bool isdna,istargetdna;//keep track of the previous settings used by oligowalk in the registry
	int option,length,usesub;
	double c;
	bool singleinsert,savefile,dynsavefile;

	CMainFrame* pMainFrame;
	bool clockwise;
	void FoldRNA(char *filename);
	void FoldDNA(char *filename);
	void BoxPlot(CString filename);
	void DotPlot(CString filename);
	void openpfwindow(CFoldDoc *pFoldDocument);
	void DisplayColorKeyWindow();
	void DisplaySHAPEKeyWindow();
	void openbimolpfwindow(CFoldDoc *pFoldDocument);
	void DynDotPlot(CString filename);

	
	CSplash *splash;


	//declarations for command line:


	int i,**table;
	char sequencefilename[maxfil],ctfilename[maxfil],sequencefilename2[maxfil],reportfilename[maxfil],confilename[maxfil];
	datatable data,*ddata;
	bool seq1,bimol,suboptimal,constrained;
	int number,percent,window,mode,conc;
	
	
	char loop[maxfil],stackf[maxfil], tstackh[maxfil], tstacki[maxfil],
		tloop[maxfil], miscloop[maxfil], danglef[maxfil], int22[maxfil],
		int21[maxfil],coax[maxfil], tstackcoax[maxfil],coaxstack[maxfil],
		tstack[maxfil], tstackm[maxfil],triloop[maxfil],int11[maxfil],loop2[maxfil],hexaloop[maxfil],
		tstacki23[maxfil], tstacki1n[maxfil];
	rddata *hybriddata;

	/*int i,**table;
	char sequencefilename[maxfil],ctfilename[maxfil],sequencefilename2[maxfil],reportfilename[maxfil];
	datatable data,*ddata;
	bool isdna,seq1,bimol,filt,suboptimal;
	int number,percent,window,fnumber,fpercent,fwindow,length,mode,conc;
	structure ct,ct2,ct3;
	double c;
	char loop[maxfil],stackf[maxfil], tstackh[maxfil], tstacki[maxfil],
		tloop[maxfil], miscloop[maxfil], danglef[maxfil], int22[maxfil],
		int21[maxfil],coax[maxfil], tstackcoax[maxfil],coaxstack[maxfil],
		tstack[maxfil], tstackm[maxfil],triloop[maxfil],int11[maxfil],loop2[maxfil];
	rddata *hybriddata;*/








private:
	CMultiDocTemplate *pFoldTemplate;
	CMultiDocTemplate *pBiFoldTemplate;
	CMultiDocTemplate *pDrawTemplate;
	CMultiDocTemplate *pSequenceTemplate;
	CMultiDocTemplate *pEfn2Template;
	CMultiDocTemplate *pOligoTemplate;
	CMultiDocTemplate *pOligoViewTemplate;
	CMultiDocTemplate *pDynalignTemplate;
	CMultiDocTemplate *pMixmatchTemplate;
	CMultiDocTemplate *pReFoldTemplate;
	CMultiDocTemplate *pDynAlignRefoldTemplate;
	CMultiDocTemplate *pPFDPTemplate;
	CMultiDocTemplate *pDPTemplate;
	CMultiDocTemplate *pPFTemplate;
	CMultiDocTemplate *pAllFoldTemplate;
	CMultiDocTemplate *pColorKeyTemplate;
	CMultiDocTemplate *pSHAPEKeyTemplate;
	CMultiDocTemplate *pOligoScreenTemplate;
	CMultiDocTemplate *pOligosiRNATemplate;
	CMultiDocTemplate *pBimolPFTemplate;
	CMultiDocTemplate *pStochasticTemplate;
	CMultiDocTemplate *pBreakPseudoTemplate;
	CMultiDocTemplate *pMEATemplate;
	CMultiDocTemplate *pPKTemplate;
	CMultiDocTemplate *pMultilignTemplate;
	CMultiDocTemplate *pTurboFoldTemplate;
	void commandline();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CRNAstructureApp)
	public:
	virtual BOOL InitInstance();
	virtual int ExitInstance();
	//}}AFX_VIRTUAL

// Implementation
	//{{AFX_MSG(CRNAstructureApp)
	afx_msg void OnAppAbout();
	afx_msg void OnFoldDNAsinglestrand();
	afx_msg void OnFoldDNAbimolecular();
	afx_msg void OnFoldRNAbimolecular();
	afx_msg void OnFileNew();
	afx_msg void OnFileEfn2();
	afx_msg void OnFileEfn2dna();
	afx_msg void OnFileOpen();
	afx_msg void OnFileOligowalk();
	afx_msg void OnFileDynalign();
	afx_msg void OnFileMixmatch();
	afx_msg void OnHelpOpen();
	afx_msg void OnFileRefoldfromsavefile();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
	afx_msg void OnDynalignRefold();
	// Create a Dynalign Dotplot
	void OnDynalignDotplot(void);
	afx_msg void OnFileDotplot();
	afx_msg void OnFileDotplotpartitionfunction();
	afx_msg void OnFilePartitionfunctionrna();
	
	afx_msg void OnFileGenerateallsuboptimalstructures();
	afx_msg void OnFileOligoscreen();
	afx_msg void OnFileOligowalkforsirna();
	afx_msg void OnFilePartitionfunctionrnabimolecular();
	afx_msg void OnFileStochastictraceback();
	afx_msg void OnFilePartitionfunctiondna();
	afx_msg void OnFileDotplotfromtextfile();
	afx_msg void OnFileBreakrnapseudoknots();
	afx_msg void OnRNAMEA();
	afx_msg void OnRnaProbknot();
	afx_msg void OnDnaProbknot();
	afx_msg void OnRNAMultilign();
public:
	afx_msg void OnDnaPartitionfunctiondnabimolecular();
public:
	afx_msg void OnDnaGeneratealldnastructures();
public:
	afx_msg void OnDnaStochasticdnasampling();
public:
	afx_msg void OnDnaDnadynalign();
public:
	afx_msg void OnDnaPredictdnameastructure();
	afx_msg void OnDnaDnamultilign();
	afx_msg void OnRnaRnaturbofold();
	afx_msg void OnDnaDnaturbofold();
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_RNASTRUCTURE_H__56CD85E4_8D8A_11D4_9F32_00C0F02A5F5D__INCLUDED_)
