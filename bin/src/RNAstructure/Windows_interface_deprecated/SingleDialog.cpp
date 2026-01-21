// SingleDialog.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "SingleDialog.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CSingleDialog dialog


CSingleDialog::CSingleDialog(structure *CT,CWnd* pParent /*=NULL*/, int Mode /*=1*/)
	: CDialog(CSingleDialog::IDD, pParent)
{
	//{{AFX_DATA_INIT(CSingleDialog)
	m_single = 0;
	//}}AFX_DATA_INIT
	ct = CT;
	mode = Mode;
	//mode = 1 -> force single stranded
	//mode = 2 -> FMN
	//mode = 3 -> chemical modification
}


void CSingleDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CSingleDialog)
	DDX_Text(pDX, IDC_NUMBER_SINGLE, m_single);
	DDV_MinMaxInt(pDX, m_single, 0, 32000);
	//}}AFX_DATA_MAP
}

BOOL CSingleDialog::OnInitDialog( ) {

	if (mode==2) {
		SetWindowText("U in GU Pair");

	}
	else if (mode==3) {
		SetWindowText("Chemically Modified");

	}


	return CDialog::OnInitDialog();

	
}

BEGIN_MESSAGE_MAP(CSingleDialog, CDialog)
	//{{AFX_MSG_MAP(CSingleDialog)
	ON_EN_CHANGE(IDC_NUMBER_SINGLE, OnChangeNumberSingle)
	ON_BN_CLICKED(IDOK2, OnOk2)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CSingleDialog message handlers

void CSingleDialog::OnChangeNumberSingle() 
{
	UpdateData(TRUE);
	if (ct->allocated) {
		if (m_single > ct->GetSequenceLength()) {
			
			AfxMessageBox( "Nucleotide is past end of sequence!", 
			MB_OK|MB_ICONHAND);
			m_single = 0;
			UpdateData(FALSE);
		}
	}
	
}

void CSingleDialog::OnOK() 
{
	
	UpdateData(TRUE);
	if (mode==2) {
		//if (ct->ngu>=maxgu) {
		//	AfxMessageBox( "Too many FMN cleaved nucleotides specified!",
		//		MB_OK|MB_ICONHAND);
		//	return;

		//}
		//else {

			//ct->gu[ct->ngu]=m_single;
			//ct->ngu++;

			ct->AddGUPair(m_single);
			
		//}


	}
	else if (mode==3) {
		//if (ct->nmod>=maxforce) {
		//	AfxMessageBox( "Too chemically modified nucleotides specified!",
		//		MB_OK|MB_ICONHAND);
		//	return;

		//}
		//else {

			//ct->nmod++;
			//ct->mod[ct->nmod]=m_single;
			ct->AddModified(m_single);
		//}

	}
	
	else if (m_single!=0) {
		//(ct->nnopair)++;
		//if (ct->nnopair>maxforce) {
		//	AfxMessageBox( "Too many single-stranded nucleotides specified!",
		///		MB_OK|MB_ICONHAND);
		//	(ct->nnopair)--;
		//	CDialog::OnOK();
			

		//}

		//ct->nopair[ct->nnopair]=m_single; 
		ct->AddSingle(m_single);
	}
	CDialog::OnOK();
}

void CSingleDialog::OnCancel() 
{
	// TODO: Add extra cleanup here
	
	CDialog::OnCancel();
}

void CSingleDialog::OnOk2() 
{
	UpdateData(TRUE);
	if (mode==2) {
		//if (ct->ngu>=maxgu) {
		//	AfxMessageBox( "Too many FMN cleaved nucleotides specified!",
		//		MB_OK|MB_ICONHAND);
		//	return;

		//}
		//else {

			//ct->gu[ct->ngu]=m_single;
			//ct->ngu++;
			ct->AddGUPair(m_single);
			
		//}


	}
	else if (mode==3) {
		//if (ct->nmod>=maxforce) {
		//	AfxMessageBox( "Too chemically modified nucleotides specified!",
		//		MB_OK|MB_ICONHAND);
		//	return;

		//}
		//else {

			//ct->nmod++;
			//ct->mod[ct->nmod]=m_single;
			ct->AddModified(m_single);
		//}

	}
	
	else if (m_single!=0) {
		//(ct->nnopair)++;
		//if (ct->nnopair>maxforce) {
		//	AfxMessageBox( "Too many single-stranded nucleotides specified!",
		//		MB_OK|MB_ICONHAND);
		//	(ct->nnopair)--;
		//	CDialog::OnOK();
			

		//}
		ct->AddSingle(m_single);
		//ct->nopair[ct->nnopair]=m_single; 
	}

	m_single = 0;
	UpdateData(FALSE);

	
}
