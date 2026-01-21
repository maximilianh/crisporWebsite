// RegionalNMR.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "RegionalNMR.h"


// CRegionalNMR dialog

IMPLEMENT_DYNAMIC(CRegionalNMR, CDialog)
CRegionalNMR::CRegionalNMR(CWnd* pParent /*=NULL*/)
	: CDialog(CRegionalNMR::IDD, pParent)
{
}

CRegionalNMR::CRegionalNMR(structure *CT, CWnd* pParent /*=NULL*/)
	: CDialog(CRegionalNMR::IDD, pParent)
	, min_g_or_u(0)
	, min_gu_pair(0)
	, neighbors(_T(""))
	, start(0)
	, stop(0)


{
	ct = CT;
}

CRegionalNMR::~CRegionalNMR()
{
}

void CRegionalNMR::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT1, min_g_or_u);
	DDX_Text(pDX, IDC_EDIT2, min_gu_pair);
	DDX_Text(pDX, IDC_EDIT3, neighbors);
	DDX_Text(pDX, IDC_START, start);
	DDX_Text(pDX, IDC_STOP, stop);
}


BEGIN_MESSAGE_MAP(CRegionalNMR, CDialog)
	ON_BN_CLICKED(IDOK, OnBnClickedOk)
	ON_BN_CLICKED(IDCANCEL, OnBnClickedCancel)
END_MESSAGE_MAP()


// CRegionalNMR message handlers


void CRegionalNMR::OnBnClickedOk()
{

	int i;
	char *string;
	//process the provided data
	UpdateData(TRUE);
	ct->rmin_gu[ct->nregion] = min_gu_pair;
	ct->rmin_g_or_u[ct->nregion] = min_g_or_u;
	ct->start[ct->nregion] = start;
	ct->stop[ct->nregion] = stop;
	
	
	//parse the string of neighbors
	string = new char[neighbors.GetLength()+1];
	//strcpy(string,);
	string = strtok(neighbors.GetBuffer(1),",");
	

	while (string) {
		//encode the data
		for (i=0;((unsigned int) i)<strlen(string)&&i<maxneighborlength;i++) {
			if (string[i]=='A'||string[i]=='a') ct->rneighbors[ct->nregion][ct->nneighbors][i]=1;
			else if (string[i]=='C'||string[i]=='c') ct->rneighbors[ct->nregion][ct->nneighbors][i]=2;
			else if (string[i]=='G'||string[i]=='g') ct->rneighbors[ct->nregion][ct->nneighbors][i]=3;
			else if (string[i]=='U'||string[i]=='u'||string[i]=='T'||string[i]=='t') ct->rneighbors[ct->nregion][ct->nneighbors][i]=4;
		}
		ct->rneighbors[ct->nregion][ct->nneighbors][strlen(string)] = 0;//zero indicates last neighbor was reached
		
		(ct->rnneighbors[ct->nregion])++;
		string = strtok(NULL,",");
	}
	ct->nregion++;
	OnOK();
}

void CRegionalNMR::OnBnClickedCancel()
{
	// TODO: Add your control notification handler code here
	OnCancel();
}

