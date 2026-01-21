// MultilignView.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "MultilignView.h"
#include "common_windows_actions.h"
#include "temp_Dialog.h"

//This is the backend functions, set up to run in a seperate thread
UINT MultilignProc( LPVOID pParam ){   
	CMultilignObject* pObject = (CMultilignObject*) pParam;
	MultilignView* pView = (MultilignView*) pObject->parent;

	Multilign_object *multilign;
	int error;


	//change to the datapath directory
	//_chdir(pObject->datapath);

	//Change the DATAPATH environment variable:
	setDataPath(pObject->datapath); //This is needed for the phmm

	//allocate the instance of multilign
	multilign = new Multilign_object(*pObject->filenames,pObject->ISRNA,pObject->progress);
	error = multilign->GetErrorCode();
	if (error!=0) {
		
		pObject->errormessage = multilign->GetErrorMessage(error).c_str();
		delete multilign;

		::PostMessage(pView->m_hWnd,ID_MAERROR,0,0);
	}

	multilign->SetMaxPairs(pObject->maxpairs);
	multilign->SetIterations(pObject->iterations);
	multilign->SetMaxDsv(pObject->maxdsvchange);
	multilign->SetTemperature(pObject->Temperature);

	error = multilign->ProgressiveMultilign(
          1,
          true, true,
          pObject->maxstruct, 
          pObject->structurewindow, pObject->alignmentwindow, 
          pObject->maxpercent, 
          -99, 
          pObject->gap, 
          pObject->singlebp, 
          30, 
          false);
	//1 = 1 processor, -99 = M, true = generate save, true = write alignments, 30 = single fold suboptimal %, final false = non-local 

	if (error!=0) {
		pObject->errormessage = multilign->GetErrorMessage(error).c_str();
		delete multilign;

		::PostMessage(pView->m_hWnd,ID_MAERROR,0,0);
		
	}

	multilign->WriteAlignment(pObject->alignment);

	//If !generatesave, then cleanup all the intermediate files
	//NOTE:  This does not work because of file locking.
	//if (!pObject->generatesave) {
	//	multilign->CleanupIntermediateFiles();
	//}


	//clean up
	delete multilign;
	
	//Send the completion message
	::PostMessage(pView->m_hWnd,ID_MULTILIGNDONE,0,0);
	
	return 0;   // thread completed successfully

}

// MultilignView

IMPLEMENT_DYNCREATE(MultilignView, CFormView)

MultilignView::MultilignView()
	: CFormView(MultilignView::IDD)
	, m_sequence(_T(""))
	, m_ct(_T(""))
	, m_alignment(_T(""))
	, m_maxpercent(20)
	, m_maxstructures(20)
	, m_structurewindow(0)
	, m_alignmentwindow(0)
	, m_gap(0.4)
	, m_singlebp(TRUE)
	, m_generatesave(TRUE)
	, m_iterations(2)
	, m_maxpairs(0)
	, m_maxdsvchange(1)
	, m_deletenumber(1)
	, Temperature(310.15)
	, started (false)
{

}

MultilignView::~MultilignView()
{
}

CMultilignDoc *MultilignView::GetMultilignDocument() {
	
	return ((CMultilignDoc*) GetDocument());	

}

void MultilignView::OnInitialUpdate() {

	
	//if (*GetMultilignDocument()->checksave) {
	//	savefilecheck = TRUE;
	//}
	//else savefilecheck = FALSE;
	

	ResizeParentToFit();
	UpdateData(FALSE);

}

void MultilignView::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_SEQUENCENAME, m_sequence);
	DDX_Text(pDX, IDC_CTNAME, m_ct);
	DDX_Text(pDX, IDC_ALIGNNAME, m_alignment);
	DDX_Text(pDX, IDC_PERCENT2, m_maxpercent);
	DDV_MinMaxDouble(pDX, m_maxpercent, 0, 100.0);
	DDX_Text(pDX, IDC_NUMBER2, m_maxstructures);
	DDX_Text(pDX, IDC_WINDOW2, m_structurewindow);
	DDX_Text(pDX, IDC_ALIGNWINDOW, m_alignmentwindow);
	DDX_Text(pDX, IDC_DGGAP, m_gap);
	DDX_Check(pDX, IDC_SINGLEINSERT, m_singlebp);
	DDX_Check(pDX, IDC_SAVECHECK, m_generatesave);
	DDX_Text(pDX, IDC_CTNAME2, m_iterations);
	DDV_MinMaxInt(pDX, m_iterations, 1, 32000000);
	DDX_Text(pDX, IDC_ALIGNMENTNAME, m_maxpairs);
	DDX_Text(pDX, IDC_ALIGNMENTNAME2, m_maxdsvchange);
	DDX_Text(pDX, IDC_DELETENUMBER, m_deletenumber);
	DDX_Control(pDX, IDC_SEQLIST, C_SequenceList);
}

BEGIN_MESSAGE_MAP(MultilignView, CFormView)
	ON_MESSAGE(ID_MULTILIGNDONE, &MultilignView::Done)
	ON_MESSAGE(ID_ERRORREADINGDATA, &MultilignView::DataError)
	ON_MESSAGE(ID_MAERROR, &MultilignView::Error)
	ON_BN_CLICKED(IDC_SEQUENCEBUTTON, &MultilignView::OnBnClickedSequencebutton)
	ON_BN_CLICKED(IDC_CTBUTTON, &MultilignView::OnBnClickedCtbutton)
	ON_BN_CLICKED(IDC_ADD2, &MultilignView::OnBnClickedAdd2)
	ON_BN_CLICKED(IDC_ALIGNMENTBUTTON, &MultilignView::OnBnClickedAlignmentbutton)
	ON_BN_CLICKED(IDC_START, &MultilignView::OnBnClickedStart)
	ON_BN_CLICKED(IDC_DELETE, &MultilignView::OnBnClickedDelete)
	ON_COMMAND(ID_TEMPERATURE, &MultilignView::OnTemperature)
END_MESSAGE_MAP()


// MultilignView diagnostics

#ifdef _DEBUG
void MultilignView::AssertValid() const
{
	CFormView::AssertValid();
}

#ifndef _WIN32_WCE
void MultilignView::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}
#endif
#endif //_DEBUG





// MultilignView message handlers


//The user is selecting a sequence and needs a filedialog
void MultilignView::OnBnClickedSequencebutton()
{
	if (GetSequenceDialog(&m_sequence, GetMultilignDocument()->startpath)) {
			
		
		GetCTName(&m_sequence, &m_ct);

		//If this is the first sequence, also give a name to the alignment
		if (m_alignment=="") {

			char *aliname;
			short int i;
			i = m_sequence.GetLength();
			
			aliname = new char[i+4];//allocate enough space so that 
															//three characters can be added 
															//to the name if necessary
			strcpy(aliname,m_sequence.GetBuffer(10));
			//count the characters to the .
			
			while(i>=0){
				
				if (aliname[i]=='.') break;
				i--;
			}
			if (i==0) i = m_sequence.GetLength();
			strcpy(aliname+i+1,"ali\0");
			m_alignment=aliname;
			
			delete[] aliname;

		}
		
		UpdateData(FALSE);
		

		
	}

}


//The user is manually choosing a ct filename
void MultilignView::OnBnClickedCtbutton()
{
	if(GetCTDialog(&m_ct, GetMultilignDocument()->startpath)) {
		UpdateData(FALSE);
		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(GetMultilignDocument()->startpath,_MAX_PATH);

	}

}


//The user is adding the current sequence and ct to the list 
void MultilignView::OnBnClickedAdd2()
{
	vector<string> oneset;


	//Make sure a sequence was specified
	if (m_sequence=="") {
		AfxMessageBox( "Please choose a sequence by clicking the \"Sequence File\"Button.", 
			MB_OK|MB_ICONHAND);
		return;

	}
	
	oneset.push_back(m_sequence.GetBuffer(0));
	oneset.push_back(m_ct.GetBuffer(0));
	oneset.push_back("");
	oneset.push_back("");

	GetMultilignDocument()->filenames.push_back(oneset);

	m_sequence = "";
	m_ct = "";

	UpdateSequenceList();

}


//Update the list of sequences displayed to the user
//This needs to be performed when a sequences is added or deleted
void MultilignView::UpdateSequenceList() {

	int i;
	CString text;
	char number[5];
	

	//Reformat the text in C_SequenceList:
	text = "";
	m_maxpairs=0;
	for (i=0;i<GetMultilignDocument()->filenames.size();++i) {
		sprintf(number,"%i",(i+1));
		text+=number;
		text+=". ";
		text+=GetMultilignDocument()->filenames[i][0].c_str();
		text+="  ";
		text+=GetMultilignDocument()->filenames[i][1].c_str();
		text+="\r\n";

		//Also calculate the mean length of sequences, which is needed for setting MaxPairs
		structure *ct;
		ct = new structure();

		ct->openseq(GetMultilignDocument()->filenames[i][0].c_str());

		//Accumulate a mean:
		m_maxpairs += ct->GetSequenceLength();

		//Set the Windows using the sie of the first sequence:
		if (i==0&&GetMultilignDocument()->filenames.size()==1) {
			
			

			if ((ct->GetSequenceLength())>1200) {
   					m_structurewindow=20;
					m_alignmentwindow=3;
			}
			else if ((ct->GetSequenceLength())>800) {
   					m_structurewindow=15;
					m_alignmentwindow=3; 
			}
			else if ((ct->GetSequenceLength())>500) {
   					m_structurewindow=11;
					m_alignmentwindow=3; 
			}
			else if ((ct->GetSequenceLength())>300) {
   					m_structurewindow=7;
					m_alignmentwindow=2; 
			}
			else if ((ct->GetSequenceLength())>120) {
   					m_structurewindow=5;
					m_alignmentwindow=1;
		            
			}
			else if ((ct->GetSequenceLength())>50) {
   					m_structurewindow=3;
					m_alignmentwindow=1;
		            
			}
			else {
				m_structurewindow=2;
				m_alignmentwindow=0;
			}
			
		}


		delete ct;


	}

	if (GetMultilignDocument()->filenames.size()>0) {
		m_maxpairs = m_maxpairs/GetMultilignDocument()->filenames.size();
	}

	C_SequenceList.SetWindowText(text);



	UpdateData(FALSE);
}


//User is overriding the automatic choice of alignment filename
void MultilignView::OnBnClickedAlignmentbutton()
{
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".ct",NULL,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"Alignment Files (*.ali)|*.ali|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		m_alignment = (filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);		

	}
	delete filedialog;
}

//The user has started the calculation
void MultilignView::OnBnClickedStart()
{
	
	//if the calculation has already been started, don't start it again!
	if (started) return;
	started = true;
	
	//Check that the sequences were chosen:
	if(GetMultilignDocument()->filenames.size()<3) {
		AfxMessageBox( "Please choose at least three sequences for this calculation.", 
			MB_OK|MB_ICONHAND);
		return;

	}

	//Get the latest values entered by the user
	UpdateData(TRUE);


	//TProgressDialog *progress;
	progress = new TProgressDialog();//parent,"Folding the RNA..."
	CRect *rect;
	rect = new CRect(10,40,210,120);
	
	if (GetMultilignDocument()->ISRNA)
		progress->Create(NULL,"Folding the RNA...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);//(IDD_PROGRESSDIALOG);//
	else progress->Create(NULL,"Folding the DNA...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);
	delete rect;

	
	object = new CMultilignObject;
	object->progress = progress;
    object->parent = this;
	object->filenames=&GetMultilignDocument()->filenames;
	object->alignment=m_alignment.GetBuffer(10);
	object->maxpercent=m_maxpercent;
	object->maxstruct=m_maxstructures;
	object->structurewindow=m_structurewindow;
	object->alignmentwindow=m_alignmentwindow;
	object->gap=m_gap;
	object->singlebp=m_singlebp;
	object->generatesave=m_generatesave;
	object->iterations=m_iterations;
	object->maxpairs=m_maxpairs;
	object->maxdsvchange=m_maxdsvchange;
	object->Temperature=Temperature;
	object->datapath=GetMultilignDocument()->datapath;
	object->ISRNA=GetMultilignDocument()->ISRNA;
	   
	AfxBeginThread(MultilignProc,object);


}


//The user is selecting to delete a sequence
void MultilignView::OnBnClickedDelete()
{
	//Update the variables
	UpdateData(TRUE);


	//Check that the number to delete makes sense:
	if (m_deletenumber > 0 && m_deletenumber <= GetMultilignDocument()->filenames.size()) {

		
		GetMultilignDocument()->filenames.erase(GetMultilignDocument()->filenames.begin() + (m_deletenumber-1));
		UpdateSequenceList();

	}
	else {
		//The number is too large or too small
		
		AfxMessageBox( "Please choose a sequence number to delete that corresponds to one of the above sequences.", 
			MB_OK|MB_ICONHAND);
		return;

		

	}


}


//The seperate thread is done and sent the message back.
LRESULT MultilignView::Done(WPARAM wParam, LPARAM lParam) {
	CDialog* finished;
	int i;


	
	//offer to display the predicted structures
	finished = new CDialog(IDD_FINISHED,this);

	if(finished->DoModal()==IDOK) {
		//draw all the structures	
		for (i=0;i<GetMultilignDocument()->filenames.size();++i) {
			((CRNAstructureApp*) GetMultilignDocument()->pMainFrame)->Draw(GetMultilignDocument()->filenames[i][1].c_str());
		}

	}
	delete finished;

	delete object;
	
	progress->SendMessage(WM_CLOSE);
	
	GetMultilignDocument()->Frame->SendMessage(WM_CLOSE);
	
	return 0;


}


//User has chosen to change the temperature
void MultilignView::OnTemperature()
{
	
	//Allow the user to specify a new temperature.
	CTemp_Dialog *temp;
	temp=new CTemp_Dialog(&Temperature);

	temp->DoModal();
	delete temp;


}

//The seperate thread is sending back the message that it had an error reading the datatables.
LRESULT MultilignView::DataError(WPARAM wParam, LPARAM lParam) {
	
	AfxMessageBox( "A thermodynamic data file could not be found!\n\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
			MB_OK|MB_ICONHAND);
	
	delete object;
	progress->SendMessage(WM_CLOSE);
	
	GetMultilignDocument()->Frame->SendMessage(WM_CLOSE);
	
	return 0;


}


//The seperate thread is sending back the message that it had an error reading the datatables.
LRESULT MultilignView::Error(WPARAM wParam, LPARAM lParam) {
	


	AfxMessageBox( object->errormessage, 
			MB_OK|MB_ICONHAND);
	
	delete object;
	progress->SendMessage(WM_CLOSE);
	
	GetMultilignDocument()->Frame->SendMessage(WM_CLOSE);
	
	return 0;


}