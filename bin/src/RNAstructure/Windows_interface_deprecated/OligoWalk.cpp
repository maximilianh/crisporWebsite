// OligoWalk.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "OligoWalk.h"
#include "../src/structure.h"
#include "../src/algorithm.h"
#include "temp_Dialog.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#include <direct.h>






//This is the backend functions for OligoWalk, set up to run in a seperate thread
UINT OligoProc( LPVOID pParam ){    
	COligoObject* oligoobject = (COligoObject*) pParam;
	COligoWalk* app = (COligoWalk*) oligoobject->parent; 
	siPREFILTER *prefilter;
	int test=-1;

	//use a character array to turn off SHAPE:
	char shapefile='\0';
    

	//Although prefiltering isn't used in GUI, it needs to be initialized and passed to function.
	//Note that the true indicates that empirical scores would be used a
	prefilter = new siPREFILTER(oligoobject->data,oligoobject->dhdata,oligoobject->useprefilter,true,oligoobject->ct.GetSequenceLength() - oligoobject->length + 2,oligoobject->isdna);
	

	//note that foldsize is temporarilly set to zero so that there is no folding size limit
	//note that distance is set to zero so that there is no maximum pairping distance
	//note that test is set to -1 so there is no testing
	//note that write is set to FALSE to tunr off writing
	olig(oligoobject->isdna, oligoobject->option, &oligoobject->ct, oligoobject->length, 
		oligoobject->c, oligoobject->table, oligoobject->numofsubstructures, oligoobject->data, *oligoobject->ddata, 
		oligoobject->hybriddata,oligoobject->usesub,oligoobject->update,oligoobject->helixstack,
		oligoobject->start,oligoobject->stop,prefilter, 0, 0, &shapefile,&test,false);

	if (oligoobject->siRNA) {
		filterbysirna(&oligoobject->ct,oligoobject->table,oligoobject->length,&oligoobject->data,oligoobject->mask,oligoobject->asuf,oligoobject->tofe,oligoobject->fnnfe);
	}


	//note that foldsize is temporarilly set to zero here.
	report(oligoobject->outputfile, &oligoobject->ct, oligoobject->table, oligoobject->numofsubstructures, 
		oligoobject->length, oligoobject->isdna, oligoobject->c, oligoobject->usesub,oligoobject->start,oligoobject->stop,prefilter,0,
		oligoobject->mask,oligoobject->asuf,oligoobject->tofe,oligoobject->fnnfe,false);




	delete prefilter;

	
	::PostMessage(app->m_hWnd,ID_OLIGODONE,0,0);	

	

	return 0;   // thread completed successfully
}


/////////////////////////////////////////////////////////////////////////////
// COligoWalk dialog


IMPLEMENT_DYNCREATE(COligoWalk, CFormView)

COligoWalk::COligoWalk(UINT i)
	: CFormView(i)
{

	

	ctloaded = false;
	

	//{{AFX_DATA_INIT(COligoWalk)
	m_length = 0;
	m_start = 0;
	m_stop = 0;
	m_ctname = _T("");
	m_reportname = _T("");
	m_concentration = 0;
	//}}AFX_DATA_INIT
	started =false;
}


void COligoWalk::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
#include "commonoligowalkdataexchange.h"
}


BEGIN_MESSAGE_MAP(COligoWalk, CFormView)
#include "commonoligowalkmessagemap.h"
	
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// COligoWalk message handlers

void COligoWalk::OnInitialUpdate() {
	
	CComboBox *m_units;
	CButton *button;
	
	m_length = (GetOligoDocument()->oligoobject->length);

	//Create the Spin Control
	CWnd* pWnd = GetDlgItem( IDC_UPDOWNLENGTH );
	CRect rect;
	pWnd->GetWindowRect( &rect );
	ScreenToClient( &rect );

	m_LengthSpin.Create( WS_VISIBLE|WS_CHILD/*|dwStyles*/, rect, this, IDC_LENGTHSPIN );
	m_LengthSpin.SetRange( 1, 200 );  // Sends UDM_SETRANGE
	m_LengthSpin.SetPos(m_length);


	pWnd = GetDlgItem( IDC_UPDOWNSTART );
	
	pWnd->GetWindowRect( &rect );
	ScreenToClient( &rect );

	m_StartSpin.Create( WS_VISIBLE|WS_CHILD/*|dwStyles*/, rect, this, IDC_STARTSPIN );
	m_StartSpin.SetRange( 1, 200 );  // Sends UDM_SETRANGE
	m_StartSpin.SetPos(1);

	pWnd = GetDlgItem( IDC_UPDOWNSTOP );
	
	pWnd->GetWindowRect( &rect );
	ScreenToClient( &rect );

	m_StopSpin.Create( WS_VISIBLE|WS_CHILD/*|dwStyles*/, rect, this, IDC_STOPSPIN );
	m_StopSpin.SetRange( 1, 200 );  // Sends UDM_SETRANGE
	m_StopSpin.SetPos(1);

	m_units = (CComboBox*) GetDlgItem (IDC_CUNITS);
	m_units->AddString("mM");
	m_units->AddString("uM");
	m_units->AddString("nM");
	m_units->AddString("pM");


	if ((GetOligoDocument()->oligoobject->c)<=1e-12) {
		m_concentration = (GetOligoDocument()->oligoobject->c)/1e-12;
		m_units->SetWindowText("pM");
	}
	else if ((GetOligoDocument()->oligoobject->c)<=1e-9) {
		m_concentration = (GetOligoDocument()->oligoobject->c)/1e-9;
		m_units->SetWindowText("nM");
	}
	else if ((GetOligoDocument()->oligoobject->c)<=1e-6) {
		m_concentration = (GetOligoDocument()->oligoobject->c)/1e-6;
		m_units->SetWindowText("uM");
	}
	else {
		m_concentration = (GetOligoDocument()->oligoobject->c)/1e-3;
		m_units->SetWindowText("mM");

	}



	if ((GetOligoDocument()->oligoobject->option)==1) {
		button = (CButton*) GetDlgItem( IDC_OPTION1 );
		button->SetCheck(1);
	}
	else if ((GetOligoDocument()->oligoobject->option)==2) {
		button = (CButton*) GetDlgItem( IDC_OPTION2 );
		button->SetCheck(1);
	}
	else {
		button = (CButton*) GetDlgItem( IDC_OPTION3 );
		button->SetCheck(1);
	}
	

	//Note that the DNA/RNA checks are invisible in OligoWalk for siRNA --
	//RNA is assumed
	if ((GetOligoDocument()->oligoobject->isdna)) {
		button = (CButton*) GetDlgItem( IDC_DNA );
		button->SetCheck(1);

	}
	else {
		button = (CButton*) GetDlgItem( IDC_RNA );
		button->SetCheck(1);

	}

	
	UpdateData(FALSE);
	ResizeParentToFit();
	
	return;
}

void COligoWalk::OnCtbutton() 
{
	CFileDialog *filedialog;
	int i;
	char *report;
	
	

	
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"CT Files (*.ct)|*.ct|Sequence Files (*.seq)|*.seq||");

	
	filedialog->m_ofn.lpstrInitialDir=GetOligoDocument()->oligoobject->startpath;
	
	if (filedialog->DoModal()==IDOK) {
		
		m_ctname=(filedialog->GetPathName()).GetBuffer(30);
		i = m_ctname.GetLength();
		
		report = new char[i+4];//allocate enough space so that 
														//three characters can be added 
														//to the name if necessary
		strcpy(report,m_ctname.GetBuffer(10));
		//count the characters to the .
		
		while(i>=0){
			
			if (report[i]=='.') break;
			i--;
		}
		if (i==0) i = m_ctname.GetLength();
		strcpy(report+i+1,"rep\0");
		m_reportname=report;
		
		delete[] report;//fix this?

		GetOligoDocument()->oligoobject->ct.openct(m_ctname.GetBuffer(0));

		//for purposes of putting a nice title in the OligoView Window:
		//strncpy(GetOligoDocument()->oligoobject->ct.ctlabel[1],m_ctname.GetBuffer(0),ctheaderlength-1);
		GetOligoDocument()->oligoobject->ct.SetCtLabel(m_ctname.GetBuffer(0),1);

		m_start = 1;
		m_stop = GetOligoDocument()->oligoobject->ct.GetSequenceLength();
		ctloaded = true;

		if (GetOligoDocument()->oligoobject->ct.GetSequenceLength()<m_length) 
			m_length = GetOligoDocument()->oligoobject->ct.GetSequenceLength();

		UpdateData(FALSE);
		
		//now store the path in Startpath so that the program can start here next time:
		//_getcwd((GetOligoDocument()->oligoobject->startpath),_MAX_PATH);

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
		strncpy(GetOligoDocument()->oligoobject->startpath,path.GetBuffer(1),i);
		*(GetOligoDocument()->oligoobject->startpath + i) ='\0';

		
		

	}
	delete filedialog;

	
}

void COligoWalk::OnReportbutton() 
{
	// User is explicitly naming the report file:
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".rep",NULL,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"Report Files (*.rep)|*.rep|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		m_reportname=(filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);		

	}
	delete filedialog;
	
}

void COligoWalk::OnLengthSpin(NMHDR* pNMHDR, LRESULT* pResult) 
{
	NM_UPDOWN* pNMUpDown = (NM_UPDOWN*)pNMHDR;
	
	
	
	if(ctloaded) {
		m_length = m_length + pNMUpDown->iDelta;
	
		if (m_length>=GetOligoDocument()->oligoobject->ct.GetSequenceLength()) { 
			m_length=GetOligoDocument()->oligoobject->ct.GetSequenceLength();
		

		}
	
		UpdateData(FALSE);
	}

	*pResult = 0;
}

void COligoWalk::OnStartSpin(NMHDR* pNMHDR, LRESULT* pResult) 
{
	NM_UPDOWN* pNMUpDown = (NM_UPDOWN*)pNMHDR;
	
	
	if (ctloaded) {
	
		m_start = m_start+ pNMUpDown->iDelta;
	
		if (m_start>=GetOligoDocument()->oligoobject->ct.GetSequenceLength()) { 
			m_start=GetOligoDocument()->oligoobject->ct.GetSequenceLength();
		

		}
		else if (m_start<1) {
			m_start = 1;

		}
	
		UpdateData(FALSE);
	}

	*pResult = 0;
}

void COligoWalk::OnStopSpin(NMHDR* pNMHDR, LRESULT* pResult) 
{
	NM_UPDOWN* pNMUpDown = (NM_UPDOWN*)pNMHDR;
	
	
	
	
	if (ctloaded) {
	
		m_stop = m_stop+ pNMUpDown->iDelta;
	
		if (m_stop>=GetOligoDocument()->oligoobject->ct.GetSequenceLength()) { 
			m_stop=GetOligoDocument()->oligoobject->ct.GetSequenceLength();
		

		}
		else if (m_stop<1) {
			m_stop = 1;

		}
	
		UpdateData(FALSE);
	}

	*pResult = 0;
}

void COligoWalk::OnStart() 
{
	CButton *button,*button2;
	UpdateData(TRUE);
	CComboBox *m_units;
	CString temp;


	//a trap to make sure the calculation isn't started more than once
	if (started) return;
	else started = true;

	if (m_start>m_stop) {
		MessageBox("Please change the limits of the walk so that stop is larger than start.");
		return;
	}

	if (!ctloaded) {
		MessageBox("Please specify a CT file name.");
		return;
	}

	button = (CButton*) GetDlgItem(IDC_DNA);
	if (button->GetCheck()) (GetOligoDocument()->oligoobject->isdna) = true;
	else (GetOligoDocument()->oligoobject->isdna) = false;


	//DNA targets are not yet supported
	(GetOligoDocument()->oligoobject->istargetdna) = false;

	
	//option 1 - break local structure
	//option 2 - refold whole RNA
	//option 3 - no local structure considered
	button = (CButton*) GetDlgItem(IDC_OPTION1);
	button2 = (CButton*) GetDlgItem(IDC_OPTION2);
	if (button->GetCheck()) (GetOligoDocument()->oligoobject->option) = 1;
	else if (button2->GetCheck()) (GetOligoDocument()->oligoobject->option) = 2;
	else (GetOligoDocument()->oligoobject->option) = 3;
	
	strcpy(GetOligoDocument()->oligoobject->outputfile,m_reportname.GetBuffer(0));
	
	(GetOligoDocument()->oligoobject->length)= m_length;

	button = (CButton*) GetDlgItem(IDC_SUBOPTIMAL);
	if (button->GetCheck()) (GetOligoDocument()->oligoobject->usesub) = 3;///!This option needed to be changed to work with John Lu's code.
	else (GetOligoDocument()->oligoobject->usesub) = 0;
	
	(GetOligoDocument()->oligoobject->start) = m_start;
	(GetOligoDocument()->oligoobject->stop) = m_stop;


	//c = C;
	(GetOligoDocument()->oligoobject->c) = m_concentration;
	m_units = (CComboBox*) GetDlgItem (IDC_CUNITS);
	m_units->GetWindowText(temp);
	if (temp == "uM") (GetOligoDocument()->oligoobject->c) = ((GetOligoDocument()->oligoobject->c))/1000000;
	else if (temp == "nM") (GetOligoDocument()->oligoobject->c) = ((GetOligoDocument()->oligoobject->c))/1000000000;
	else if (temp == "mM") (GetOligoDocument()->oligoobject->c) = ((GetOligoDocument()->oligoobject->c))/1000;
	else if (temp == "pM") (GetOligoDocument()->oligoobject->c) = ((GetOligoDocument()->oligoobject->c))/1000000000000;

	if (GetOligoDocument()->oligoobject->siRNA) {
		//This is siRNA, so fetch those parameters
		GetOligoDocument()->oligoobject->asuf = m_asuf;
		GetOligoDocument()->oligoobject->tofe = m_tofe;
		GetOligoDocument()->oligoobject->fnnfe = m_fnnfe;

	}
	else {
		GetOligoDocument()->oligoobject->asuf = 0;
		GetOligoDocument()->oligoobject->tofe = 0;
		GetOligoDocument()->oligoobject->fnnfe = 0;

	}


	//For now, do not use the prefilter (for siRNA design)
	GetOligoDocument()->oligoobject->useprefilter=0;

	
	if(!GetOligoDocument()->oligoobject->AllocateTable()) {
		
		
		AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
			delete GetOligoDocument()->oligoobject;
			GetOligoDocument()->Frame->SendMessage(WM_CLOSE);
			return;

	}

	else {

		GetOligoDocument()->oligoobject->update = new TProgressDialog();//parent,"Folding the RNA..."
		CRect *rect;
		rect = new CRect(10,40,210,120);
	
		
		GetOligoDocument()->oligoobject->update->Create(NULL,"Calculating OligoWalk...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);//(IDD_PROGRESSDIALOG);//
		
		delete rect;

		CFormView* conversion = (CFormView *) this;

		GetOligoDocument()->oligoobject->parent = conversion;

		      
		
		
		AfxBeginThread(OligoProc,GetOligoDocument()->oligoobject);

	}
	
	
	
	
	
	
	
	

}

void COligoWalk::OnChangeLength() 
{

	UpdateData(TRUE);
	if (ctloaded) {
		if (m_length>GetOligoDocument()->oligoobject->ct.GetSequenceLength()) {
			m_length = GetOligoDocument()->oligoobject->ct.GetSequenceLength();
			UpdateData(FALSE);
		}

	}
		
}


void COligoWalk::OnChangeStart() 
{
	UpdateData(TRUE);
	if (ctloaded) {
		if (m_start>GetOligoDocument()->oligoobject->ct.GetSequenceLength()) {
			m_start = GetOligoDocument()->oligoobject->ct.GetSequenceLength();
			UpdateData(FALSE);
		}

		if (m_start<1) {
			m_start = 1;
			UpdateData(FALSE);
		}
	}
	
}

void COligoWalk::OnChangeStop() 
{
	UpdateData(TRUE);
	if (ctloaded) {
		if (m_stop>GetOligoDocument()->oligoobject->ct.GetSequenceLength()) {
			m_stop = GetOligoDocument()->oligoobject->ct.GetSequenceLength();
			UpdateData(FALSE);
		}

		if (m_stop<1) {
			m_stop = 1;
			UpdateData(FALSE);
		}
	}
	
}


COligoDoc* COligoWalk::GetOligoDocument() {
	
	return ((COligoDoc*) GetDocument());	

}

LRESULT COligoWalk::DisplayOligoWalk(WPARAM wParam, LPARAM lParam) {
	

	CRNAstructureApp* app = (CRNAstructureApp*) GetOligoDocument()->pParent; 
	app->OligoWalk(GetOligoDocument()->oligoobject);

	GetOligoDocument()->oligoobject->update->SendMessage (WM_CLOSE);

	GetOligoDocument()->Frame->SendMessage(WM_CLOSE);
	return 0;

}

void COligoWalk::OnTemperature()
{
	//Allow the user to specify a new temperature.
	CTemp_Dialog *temp;
	temp=new CTemp_Dialog(&(GetOligoDocument()->oligoobject->T));

	temp->DoModal();
	delete temp;
	


}