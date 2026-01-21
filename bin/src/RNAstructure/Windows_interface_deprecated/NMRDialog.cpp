// NMRDialog.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "NMRDialog.h"
#include <string.h>


// NMRDialog dialog

IMPLEMENT_DYNAMIC(NMRDialog, CDialog)
NMRDialog::NMRDialog(CWnd* pParent /*=NULL*/)
	: CDialog(NMRDialog::IDD, pParent)
	, min_g_or_u(0)
	, min_gu_pair(0)
	, neighbors(_T(""))


{
}

NMRDialog::NMRDialog(structure *CT, CWnd* pParent /*=NULL*/)
	: CDialog(NMRDialog::IDD, pParent)
	, min_g_or_u(0)
	, min_gu_pair(0)
	, neighbors(_T(""))


{
	ct = CT;
}

NMRDialog::~NMRDialog()
{
}

void NMRDialog::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT1, min_g_or_u);
	DDX_Text(pDX, IDC_EDIT2, min_gu_pair);
	DDX_Text(pDX, IDC_EDIT3, neighbors);
}


BEGIN_MESSAGE_MAP(NMRDialog, CDialog)
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
	ON_BN_CLICKED(IDCANCEL, OnBnClickedCancel)
END_MESSAGE_MAP()


// NMRDialog message handlers

void NMRDialog::OnBnClickedOk()
{
	int i;
	char *string;
	//process the provided data
	UpdateData(TRUE);
	ct->min_gu = min_gu_pair;
	ct->min_g_or_u = min_g_or_u;
	
	//parse the string of neighbors
	string = new char[neighbors.GetLength()+1];
	//strcpy(string,);
	string = strtok(neighbors.GetBuffer(1),",");
	

	while (string) {
		//encode the data
		for (i=0;(size_t) i<strlen(string)&&i<maxneighborlength;i++) {
			if (string[i]=='A'||string[i]=='a') ct->neighbors[ct->nneighbors][i]=1;
			else if (string[i]=='C'||string[i]=='c') ct->neighbors[ct->nneighbors][i]=2;
			else if (string[i]=='G'||string[i]=='g') ct->neighbors[ct->nneighbors][i]=3;
			else if (string[i]=='U'||string[i]=='u'||string[i]=='T'||string[i]=='t') ct->neighbors[ct->nneighbors][i]=4;
		}
		ct->neighbors[ct->nneighbors][strlen(string)] = 0;//zero indicates last neighbor was reached
		
		ct->nneighbors++;
		string = strtok(NULL,",");
	}

	OnOK();
}

void NMRDialog::OnBnClickedCancel()
{
	// TODO: Add your control notification handler code here
	OnCancel();
}
