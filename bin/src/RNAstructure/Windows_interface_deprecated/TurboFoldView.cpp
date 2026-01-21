// TurboFoldView.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "TurboFoldView.h"
#include "common_windows_actions.h"
#include "temp_Dialog.h"



//This is the backend functions, set up to run in a seperate thread
UINT TurboFoldProc( LPVOID pParam ){   
	CTurboFoldObject* pObject = (CTurboFoldObject*) pParam;
	CTurboFoldView* pView = (CTurboFoldView*) pObject->parent;

	TurboFold *turbofold;
	int error;
	int i;


	//Change the DATAPATH environment variable:
	_putenv_s("DATAPATH",pObject->datapath);

	//allocate the instance of turbofold
	turbofold = new TurboFold(pObject->seqfilenames,pObject->savefilenames);
	error = turbofold->GetErrorCode();
	if (error!=0) {
		
		pObject->errormessage = turbofold->GetErrorMessage(error);
		delete turbofold;

		::PostMessage(pView->m_hWnd,ID_MAERROR,0,0);
	}

	turbofold->SetProgress(*pObject->progress);
	

	error = turbofold->SetTemperature((double) pObject->Temperature); 

	if (error!=0) {
		pObject->errormessage = turbofold->GetErrorMessage(error);
		delete turbofold;

		::PostMessage(pView->m_hWnd,ID_MAERROR,0,0);
		
	}


	error = turbofold->fold(pObject->gamma,pObject->iterations); 

	if (error!=0) {
		pObject->errormessage = turbofold->GetErrorMessage(error);
		delete turbofold;

		::PostMessage(pView->m_hWnd,ID_MAERROR,0,0);
		
	}

	//Now generate and write the structures:
	for (i=0;i<pObject->seqfilenames->size();i++) {
		if (pObject->mode==1) {
			//MEA
			turbofold->MaximizeExpectedAccuracy(i+1,pObject->maxpercent,pObject->maxstruct,pObject->window,pObject->meagamma);


		}
		else if (pObject->mode==2) {
			//Threshold
			turbofold->PredictProbablePairs(i+1,pObject->probability);


		}

		else { //mode==3
			//TurboKnot/ProbKnot
			turbofold->ProbKnot(i+1, pObject->knotiterations, pObject->minlength);


		}

		turbofold->WriteCt(i+1,(*pObject->outputctfilenames)[i].c_str());

	}




	//clean up
	delete turbofold;
	
	//Send the completion message
	::PostMessage(pView->m_hWnd,ID_TURBOFOLDDONE,0,0);
	
	return 0;   // thread completed successfully

}

// CTurboFoldView

IMPLEMENT_DYNCREATE(CTurboFoldView, CFormView)

CTurboFoldView::CTurboFoldView()
	: CFormView(CTurboFoldView::IDD)
	, m_sequence(_T(""))
	, m_ct(_T(""))
	, m_mea(true)
	, m_pt(false)
	, m_gamma(0.3)
	, m_iterations(3)
	, m_percent(50)
	, m_maxstructures(1000)
	, m_window(5)
	, m_meagamma(1.0)
	, m_threshold(0)
	, m_delete(1)
	, m_pk(false)
	, m_knotiterations(1)
	, m_minlength(3)
{
	Temperature=310.15;
	started=false;
}

CTurboFoldView::~CTurboFoldView()
{
}

void CTurboFoldView::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_SEQUENCENAME, m_sequence);
	DDX_Text(pDX, IDC_CTNAME, m_ct);
	DDX_Text(pDX, IDC_GAMMA, m_gamma);
	DDX_Text(pDX, IDC_ITERATIONS, m_iterations);
	DDV_MinMaxInt(pDX, m_iterations, 1, 1024);
	DDX_Text(pDX, IDC_PERCENT2, m_percent);
	DDV_MinMaxDouble(pDX, m_percent, 0, 100);
	DDX_Text(pDX, IDC_NUMBER2, m_maxstructures);
	DDV_MinMaxInt(pDX, m_maxstructures, 1, 32000);
	DDX_Text(pDX, IDC_WINDOW2, m_window);
	DDV_MinMaxInt(pDX, m_window, 0, 32000);
	DDX_Text(pDX, IDC_MEAGAMMA, m_meagamma);
	DDX_Text(pDX, IDC_THRESHOLD, m_threshold);
	DDV_MinMaxDouble(pDX, m_threshold, 0, 100);
	DDX_Text(pDX, IDC_DELETENUMBER, m_delete);
	DDV_MinMaxInt(pDX, m_delete, 1, 32000);
	DDX_Control(pDX, IDC_SEQLIST, C_SequenceList);
	DDX_Text(pDX, IDC_TURBOKNOTITERATIONS, m_knotiterations);
	DDX_Text(pDX, IDC_ML, m_minlength);
}

BEGIN_MESSAGE_MAP(CTurboFoldView, CFormView)
	ON_MESSAGE(ID_TURBOFOLDDONE, &CTurboFoldView::Done)
	ON_MESSAGE(ID_ERRORREADINGDATA, &CTurboFoldView::DataError)
	ON_MESSAGE(ID_MAERROR, &CTurboFoldView::Error)
	ON_BN_CLICKED(IDC_SEQUENCEBUTTON, &CTurboFoldView::OnSequencebutton)
	ON_BN_CLICKED(IDC_CTBUTTON, &CTurboFoldView::OnCtbutton)
	ON_BN_CLICKED(IDC_ADD2, &CTurboFoldView::OnAdd)
	ON_BN_CLICKED(IDC_DELETE, &CTurboFoldView::OnDelete)
	ON_BN_CLICKED(IDC_START, &CTurboFoldView::OnStart)
	ON_COMMAND(ID_TEMPERATURE, &CTurboFoldView::OnTemperature)
END_MESSAGE_MAP()


// CTurboFoldView diagnostics

#ifdef _DEBUG
void CTurboFoldView::AssertValid() const
{
	CFormView::AssertValid();
}

#ifndef _WIN32_WCE
void CTurboFoldView::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}
#endif
#endif //_DEBUG


CTurboFoldDoc *CTurboFoldView::GetTurboFoldDocument() {
	
	return ((CTurboFoldDoc*) GetDocument());	

}


void CTurboFoldView::OnInitialUpdate() {

	

	ResizeParentToFit();
	UpdateData(FALSE);


	//Set the MEA button on by default
	CButton* button = (CButton*) GetDlgItem( IDC_MEA );
	button->SetCheck(1);

	

}


// CTurboFoldView message handlers


void CTurboFoldView::OnSequencebutton()
{
	if (GetSequenceDialog(&m_sequence, GetTurboFoldDocument()->startpath)) {
			
		
		GetCTName(&m_sequence, &m_ct);
	
		
		UpdateData(FALSE);
		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(GetTurboFoldDocument()->startpath,_MAX_PATH);

		
	}
}


void CTurboFoldView::OnCtbutton()
{
	if(GetCTDialog(&m_ct, GetTurboFoldDocument()->startpath)) {
		UpdateData(FALSE);
		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(GetTurboFoldDocument()->startpath,_MAX_PATH);

	}
}


void CTurboFoldView::OnAdd()
{
	
	

	//Make sure a sequence was specified
	if (m_sequence=="") {
		AfxMessageBox( "Please choose a sequence by clicking the \"Sequence File\"Button.", 
			MB_OK|MB_ICONHAND);
		return;

	}
	
	GetTurboFoldDocument()->seqfilenames.push_back(m_sequence.GetBuffer(0));
	GetTurboFoldDocument()->ctfilenames.push_back(m_ct.GetBuffer(0));


	//Now also add a pfs file to the savefilenames list:
	char *pfsname;
	short int i;
	i = m_sequence.GetLength();
	
	pfsname = new char[i+4];//allocate enough space so that 
													//three characters can be added 
													//to the name if necessary
	strcpy(pfsname,m_sequence.GetBuffer(10));
	//count the characters to the .
	
	while(i>=0){
		
		if (pfsname[i]=='.') break;
		i--;
	}
	if (i==0) i = m_sequence.GetLength();
	strcpy(pfsname+i+1,"pfs\0");
	
	

	GetTurboFoldDocument()->savefilenames.push_back(pfsname);

	delete[] pfsname;
	

	m_sequence = "";
	m_ct = "";

	UpdateSequenceList();
}

//Update the list of sequences displayed to the user
//This needs to be performed when a sequences is added or deleted
void CTurboFoldView::UpdateSequenceList() {

	int i;
	CString text;
	char number[5];
	

	//Reformat the text in C_SequenceList:
	text = "";
	
	for (i=0;i<GetTurboFoldDocument()->seqfilenames.size();++i) {
		sprintf(number,"%i",(i+1));
		text+=number;
		text+=". ";
		text+=GetTurboFoldDocument()->seqfilenames[i].c_str();
		text+="  ";
		text+=GetTurboFoldDocument()->ctfilenames[i].c_str();
		text+="\r\n";

	}
	

	C_SequenceList.SetWindowText(text);



	



	UpdateData(FALSE);
}


void CTurboFoldView::OnDelete()
{
	//Update the variables
	UpdateData(TRUE);


	//Check that the number to delete makes sense:
	if (m_delete > 0 && m_delete <= GetTurboFoldDocument()->seqfilenames.size()) {

		
		GetTurboFoldDocument()->seqfilenames.erase(GetTurboFoldDocument()->seqfilenames.begin() + (m_delete-1));
		GetTurboFoldDocument()->savefilenames.erase(GetTurboFoldDocument()->savefilenames.begin() + (m_delete-1));
		GetTurboFoldDocument()->ctfilenames.erase(GetTurboFoldDocument()->ctfilenames.begin() + (m_delete-1));
		UpdateSequenceList();

	}
	else {
		//The number is too large or too small
		
		AfxMessageBox( "Please choose a sequence number to delete that corresponds to one of the above sequences.", 
			MB_OK|MB_ICONHAND);
		return;

		

	}
}


void CTurboFoldView::OnStart()
{
	//if the calculation has already been started, don't start it again!
	if (started) return;
	started = true;
	
	//Check that the sequences were chosen:
	if(GetTurboFoldDocument()->seqfilenames.size()<2) {
		AfxMessageBox( "Please choose at least two sequences for this calculation.", 
			MB_OK|MB_ICONHAND);
		return;

	}

	//Get the latest values entered by the user
	UpdateData(TRUE);


	
	progress = new TProgressDialog();//parent,"Folding the RNA..."
	CRect *rect;
	rect = new CRect(10,40,210,120);
	
	if (GetTurboFoldDocument()->ISRNA)
		progress->Create(NULL,"Folding the RNA...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);//(IDD_PROGRESSDIALOG);//
	else progress->Create(NULL,"Folding the DNA...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);
	delete rect;


	//Set the contents of object, which will be passed to the the thread that will do the structure prediction:
	object=new CTurboFoldObject;
	object->progress=progress;
    object->parent=this;

	object->savefilenames=&GetTurboFoldDocument()->savefilenames;
	object->seqfilenames=&GetTurboFoldDocument()->seqfilenames;
	object->outputctfilenames=&GetTurboFoldDocument()->ctfilenames;

	object->Temperature=Temperature;
	
	
	if (m_mea) object->mode=1;//1 = MEA, 2= Threshold, 3 = ProbKnot
	else if (m_threshold) object->mode=2;
	else object->mode = 3;

	object->gamma=m_gamma;
	object->iterations=m_iterations;

	object->probability= m_threshold/100.0;//needs to be between 0 and 1

	object->maxpercent=m_percent;
	object->meagamma=m_meagamma;//maxperecnt is between 0 and 100
	object->maxstruct=m_maxstructures;
	object->window=m_window;

	object->datapath=GetTurboFoldDocument()->datapath;
	
	
	object->ISRNA=GetTurboFoldDocument()->ISRNA;

	object->knotiterations=m_knotiterations;
	object->minlength=m_minlength;
	
	   
	AfxBeginThread(TurboFoldProc,object);

}


//The seperate thread is done and sent the message back.
LRESULT CTurboFoldView::Done(WPARAM wParam, LPARAM lParam) {
	
	CDialog* finished;
	int i;


	
	//offer to display the predicted structures
	finished = new CDialog(IDD_FINISHED,this);

	if(finished->DoModal()==IDOK) {
		//draw all the structures	
		for (i=0;i<GetTurboFoldDocument()->ctfilenames.size();++i) {
			((CRNAstructureApp*) GetTurboFoldDocument()->pMainFrame)->Draw(GetTurboFoldDocument()->ctfilenames[i].c_str());
			((CRNAstructureApp*)GetTurboFoldDocument()->pMainFrame)->BoxPlot(GetTurboFoldDocument()->savefilenames[i].c_str());
		}

	}
	delete finished;

	delete object;
	
	progress->SendMessage(WM_CLOSE);
	
	GetTurboFoldDocument()->Frame->SendMessage(WM_CLOSE);
	
	return 0;


}


//User has chosen to change the temperature
void CTurboFoldView::OnTemperature()
{
	
	//Allow the user to specify a new temperature.
	CTemp_Dialog *temp;
	temp=new CTemp_Dialog(&Temperature);

	temp->DoModal();
	delete temp;


}

//The seperate thread is sending back the message that it had an error reading the datatables.
LRESULT CTurboFoldView::DataError(WPARAM wParam, LPARAM lParam) {
	
	AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
	
	delete object;
	
	progress->SendMessage(WM_CLOSE);
	
	GetTurboFoldDocument()->Frame->SendMessage(WM_CLOSE);
	
	return 0;


}


//The seperate thread is sending back the message that it had an error reading the datatables.
LRESULT CTurboFoldView::Error(WPARAM wParam, LPARAM lParam) {
	


	AfxMessageBox( object->errormessage, 
			MB_OK|MB_ICONHAND);
	
	delete object;
	progress->SendMessage(WM_CLOSE);
	
	GetTurboFoldDocument()->Frame->SendMessage(WM_CLOSE);
	
	return 0;


}