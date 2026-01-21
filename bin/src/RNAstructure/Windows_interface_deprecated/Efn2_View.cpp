// Efn2_View.cpp : implementation file
//

#include "stdafx.h"
#include <string>
#include <sstream>
#include "RNAstructure.h"
#include "Efn2_View.h"
#include <direct.h>
#include "../src/algorithm.h"
#include "temp_Dialog.h"
#include "../RNA_class/RNA.h"

using namespace std;


#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CEfn2_View

IMPLEMENT_DYNCREATE(CEfn2_View, CFormView)

CEfn2_View::CEfn2_View()
	: CFormView(CEfn2_View::IDD)
	, m_details(FALSE)
{
	//{{AFX_DATA_INIT(CEfn2_View)
	m_ctfilename = _T("");
	m_outfilename = _T("");
	//}}AFX_DATA_INIT
}

CEfn2_View::~CEfn2_View()
{
}

void CEfn2_View::OnInitialUpdate()
{

	//ResizeParentToFit(FALSE);
	ResizeParentToFit();
	

}

void CEfn2_View::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CEfn2_View)
	DDX_Text(pDX, IDC_CT_FILENAME, m_ctfilename);
	DDX_Text(pDX, IDC_OUT_FILENAME, m_outfilename);
	//}}AFX_DATA_MAP
	DDX_Check(pDX, IDC_DETAILS, m_details);
}


BEGIN_MESSAGE_MAP(CEfn2_View, CFormView)
	//{{AFX_MSG_MAP(CEfn2_View)
	ON_BN_CLICKED(IDC_CTFILE_EFN, OnCtfileEfn)
	ON_BN_CLICKED(IDC_OUTFILE_EFN, OnOutfileEfn)
	ON_BN_CLICKED(IDC_STARTEFN, OnStartefn)
	//}}AFX_MSG_MAP
	ON_COMMAND(ID_TEMPERATURE, OnTemperature)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CEfn2_View diagnostics

#ifdef _DEBUG
void CEfn2_View::AssertValid() const
{
	CFormView::AssertValid();
}

void CEfn2_View::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CEfn2_View message handlers

void CEfn2_View::OnCtfileEfn() 
{
		CFileDialog *filedialog;
	
	char *outname;
	short int i;

	
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"CT Files (*.ct)|*.ct||");

	
	filedialog->m_ofn.lpstrInitialDir=GetFoldDocument()->startpath;
	if (filedialog->DoModal()==IDOK) {
		//strcpy(m_sequencename.GetBuffer(10),(filedialog->GetPathName()).GetBuffer(0));
		m_ctfilename=(filedialog->GetPathName()).GetBuffer(30);
		
		
		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(GetFoldDocument()->startpath,_MAX_PATH);

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
		strncpy(GetFoldDocument()->startpath,path.GetBuffer(1),i);
		*(GetFoldDocument()->startpath + i) ='\0';



	
		i = m_ctfilename.GetLength();
		
		outname = new char[i+5];//allocate enough space so that 
														//three characters can be added 
														//to the name if necessary
		strcpy(outname,m_ctfilename.GetBuffer(10));
		//count the characters to the .
		
		while(i>=0){
			
			if (outname[i]=='.') break;
			i--;
		}
		if (i==0) i = m_ctfilename.GetLength();
		strcpy(outname+i+1,"out\0");
		m_outfilename=outname;
		
		delete[] outname;//fix this?

		UpdateData(FALSE);
		

	}
	delete filedialog;
	
}

void CEfn2_View::OnOutfileEfn() 
{
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".ct",NULL,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"Out Files (*.out)|*.out|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		m_outfilename=(filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);		

	}
	delete filedialog;
	
}

void CEfn2_View::OnStartefn() 
{
	//structure ct;
	//Switch to using RNA class (9/30/2010)
	RNA *ct;
	int i;
	string error;
	stringstream crt;
	ofstream *out;
	
	UpdateData(TRUE);	
	if (m_ctfilename==""||m_outfilename=="") {
		MessageBox("Please specify both file names.");
		return;
	}

	//Allocate the instance of RNA
	ct = new RNA(m_ctfilename.GetBuffer(1),FILE_CT,GetFoldDocument()->ISRNA);

	//Check for errors
	if(ct->GetErrorCode()!=0) {
		
		AfxMessageBox( ct->GetErrorMessage(ct->GetErrorCode()), 
				MB_OK|MB_ICONHAND);
		delete ct;
		return;

	}

	//Set the folding temperature
	ct->SetTemperature(GetFoldDocument()->T);

	//Explicitly call the read of thermodynamic parameters.  This allows the specification of the path:
	ct->ReadThermodynamic(GetFoldDocument()->Datapath);


	//Check for errors
	if(ct->GetErrorCode()!=0) {
		if(ct->GetErrorCode()==5) {
			//This is the error code when the thermodynamic parameters were not found
			AfxMessageBox( "An thermodynamic data file could not be found!\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
				MB_OK|MB_ICONHAND);
		}
		
		else {
			//There was some other error
			AfxMessageBox( ct->GetErrorMessage(ct->GetErrorCode()), 
				MB_OK|MB_ICONHAND);
		}
		delete ct;
		return;

	}

	//check to see if the temperature has been changed.
	//if (GetFoldDocument()->T<310||GetFoldDocument()->T>311) {

		//change the temperature from 310.15
	//	if (GetFoldDocument()->newtemp()==0) {
			//if newtemp returned zero, pass a warning to the user
	//		AfxMessageBox( "An enthalpy data file could not be found!\nTemperature of prediction will revert back to 37 degrees C.\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
	//			MB_OK|MB_ICONHAND);

	//	}

	//}
	//openct(&ct,m_ctfilename.GetBuffer(0));


	//Something new:  Check for pseudoknots
	for (i=1;i<=ct->GetStructureNumber();i++) {
		if (ct->ContainsPseudoknot(i)) {
			//There was a pseudoknot in structure i
			error = "Pseudoknots are not allowed.  There is a pseudoknot in structure #";
			crt << i << "\n"; 
			error+=crt.str();
			AfxMessageBox( error.c_str(), MB_OK|MB_ICONHAND);\
			delete ct;
			return;


		}

	}

	

	//process the thermodynamic details checkbox
	if (m_details) {
		//User would like the thermodynamic details:
		//efn2 (&(GetFoldDocument()->data),&ct,0,false,m_outfilename.GetBuffer(0));
		ct->WriteThermodynamicDetails(m_outfilename.GetBuffer(0));
		
	}
	else {
		//efn2 (&(GetFoldDocument()->data),&ct);
		//energyout (&ct, m_outfilename.GetBuffer(0));


		out = new ofstream(m_outfilename.GetBuffer(0));
		for (i=1;i<=ct->GetStructureNumber();i++){

			
			*out << "Structure: "<<i<<"   Energy = "<<ct->CalculateFreeEnergy(i)<<"   \n";
			//Check for errors
			if(ct->GetErrorCode()!=0) {
		
				AfxMessageBox( ct->GetErrorMessage(ct->GetErrorCode()), 
						MB_OK|MB_ICONHAND);
				delete ct;
				return;

			}
			

		}
		out->close();
		delete out;
	}

	
	MessageBox("EFN2 finished", "Calculation is complete.");
    GetFoldDocument()->Frame->SendMessage(WM_CLOSE);
	
}

CFoldDoc *CEfn2_View::GetFoldDocument() {
	
	return ((CFoldDoc*) GetDocument());	

}

void CEfn2_View::OnTemperature()
{
	//Allow the user to specify a new temperature.
	CTemp_Dialog *temp;
	temp=new CTemp_Dialog(&(GetFoldDocument()->T));

	temp->DoModal();
	delete temp;
	


}