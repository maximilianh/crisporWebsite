// Dynalign.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "Dynaligninterface.h"
#include "../src/algorithm.h"
#include "../src/dynalign.h"
#include "doubledialog.h"
#include "pairdialog.h"
#include "singledialog.h"
#include "../src/outputconstraints.h"
#include "foldview.h"
#include "alignforcedialog.h"
#include "temp_Dialog.h"
#include "../src/phmm/phmm_aln.h"
#include "../src/phmm/structure/structure_object.h"
#include "error.h"


#include <direct.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#define filter //filter out pairs from dynalign calculation
//#undef filter

using namespace std;


//This is the backend functions, set up to run in a seperate thread
UINT DynProc( LPVOID pParam ){    
	CDynObject* pObject = (CDynObject*) pParam;
	CDynalign* pView = (CDynalign*) pObject->parent; 
	short int **align;
	short i;
	bool constraints;
	bool **allowed_alignments;
	
	

	//Change the DATAPATH environment variable:
	setDataPath(pObject->datapath); //This is needed for the phmm
	
	align = new short *[pObject->maxtracebacks];//maximum number of tracebacks and next line and below at delete
	for (i=0;i<pObject->maxtracebacks;i++)  align[i] = new short [pObject->ct1->GetSequenceLength()+1];

	


	constraints = false;
	if (pObject->ct1->GetNumberofDoubles()>0) constraints=true;
	if (pObject->ct2->GetNumberofDoubles()>0) constraints=true;
	if (pObject->ct1->GetNumberofSingles()>0) constraints=true;
	if (pObject->ct2->GetNumberofSingles()>0) constraints=true;
	if (pObject->ct1->GetNumberofForbiddenPairs()>0) constraints=true;
	if (pObject->ct2->GetNumberofForbiddenPairs()>0) constraints=true;
	if (pObject->ct1->GetNumberofModified()>0) constraints=true;
	if (pObject->ct2->GetNumberofModified()>0) constraints=true;
	if (pObject->ct1->GetNumberofGU()>0) constraints=true;
	if (pObject->ct2->GetNumberofGU()>0) constraints=true;

#ifdef filter
	//This section folds the individual sequences to find insignificant pairs to be ingnored by Dynalign
	//dynamic(pObject->ct1,pObject->data,1, 0,0,NULL /*&progress*/,false, "temp1.sav");
	//dynamic(pObject->ct2,pObject->data,1, 0,0,NULL /*&progress*/,false, "temp2.sav");
	pObject->ct1->allocatetem();
	pObject->ct2->allocatetem();
	templatefromfold(pObject->ct1, pObject->data, PERCENTDOTS);
	templatefromfold(pObject->ct2, pObject->data, PERCENTDOTS);
#endif

	
	if (pObject->m == -99) {
		//allocate space in allowed_alignments, the code is using posterior probs 
			//from HMM alignment to constrain the alignment space
		allowed_alignments = new bool *[pObject->ct1->GetSequenceLength()+1];
		for (i=0;i<=pObject->ct1->GetSequenceLength();i++) {
			allowed_alignments[i] = new bool [pObject->ct2->GetSequenceLength()+1];	
	
		}

		// Needed for having nucleotide sequences as c strings.
		pObject->ct1->nucs[pObject->ct1->GetSequenceLength() + 1] = 0;
		pObject->ct2->nucs[pObject->ct2->GetSequenceLength() + 1] = 0;

		//now run HMM to determine allowed positions:
		

		//calculate_coinc_probs_env(pObject->ct1, pObject->ct2, NULL, allowed_alignments);
		//This syncs to Arif's changes
		//calculate_coinc_probs_env(pObject->ct1, pObject->ct2, NULL, allowed_alignments, pObject->datapath);
		calculate_coinc_probs_env(pObject->ct1, pObject->ct2, allowed_alignments, pObject->forcealign);

	}
	else allowed_alignments = NULL;
	if (dynalign(pObject->ct1,pObject->ct2,align,pObject->m,pObject->dggap,
		pObject->data,pObject->singleinsert,pObject->maxtracebacks,pObject->structwindow,pObject->alignwindow,
		pObject->percent,pObject->forcealign,allowed_alignments,pObject->progress,pObject->savefile,false,false,constraints)==14) {
		
			//if 14 was return, there was a traceback error that needs to be reported to the user:
		CError *error=new CError();
		error->DoModal();
		delete error;
	}
	pObject->ct1->ctout(pObject->ctoutfile);
	pObject->ct2->ctout(pObject->ctoutfile2);


	
	alignout(align,pObject->aoutfile,pObject->ct1,pObject->ct2);



	//for (i=1;i<ct1.numofbases+ct2.numofbases;i++) delete[] alignment[i];
	for (i=0;i<pObject->maxtracebacks;i++)  delete[] align[i];
	delete[] align;

	
      

	
	
	::PostMessage(pView->m_hWnd,ID_DYNDONE,0,0);
	

	return 0;   // thread completed successfully
}


/////////////////////////////////////////////////////////////////////////////
// CDynalign

IMPLEMENT_DYNCREATE(CDynalign, CFormView)

CDynalign::CDynalign()
	: CFormView(CDynalign::IDD)
	, percent(0)
	, maxtracebacks(0)
	, structwindow(0)
	, alignwindow(0)
	, savefilecheck(FALSE)
{
	//{{AFX_DATA_INIT(CDynalign)
	m_ctname = _T("");
	m_ctname2 = _T("");
	m_sequence2name = _T("");
	m_sequencename = _T("");	
	m_alignment = _T("");
	m_singleinsert = TRUE;
	m_dggap = 0.4;
	//}}AFX_DATA_INIT
	m_M = -99;//default will use posterior probability-calculated allowed alignments
	forcealign=NULL;
}

CDynalign::~CDynalign()
{

	if (forcealign!=NULL) {
		delete[] forcealign[0];
		delete[] forcealign[1];
		delete[] forcealign;

	}
}


void CDynalign::OnInitialUpdate() {

	
	if (*GetDynDocument()->checksave) {
		savefilecheck = TRUE;
		

	}
	else savefilecheck = FALSE;
	

	ResizeParentToFit();
	UpdateData(FALSE);

}

void CDynalign::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CDynalign)
	DDX_Text(pDX, IDC_CTNAME, m_ctname);
	DDX_Text(pDX, IDC_CTNAME2, m_ctname2);
	DDX_Text(pDX, IDC_SEQUENCE2NAME, m_sequence2name);
	DDX_Text(pDX, IDC_SEQUENCENAME, m_sequencename);
	DDX_Text(pDX, IDC_ALIGNMENTNAME, m_alignment);
	DDX_Check(pDX, IDC_SINGLEINSERT, m_singleinsert);
	DDX_Text(pDX, IDC_DGGAP, m_dggap);
	//}}AFX_DATA_MAP
	DDX_Text(pDX, IDC_PERCENT2, percent);
	DDV_MinMaxInt(pDX, percent, 0, 100);
	DDX_Text(pDX, IDC_NUMBER2, maxtracebacks);
	DDX_Text(pDX, IDC_WINDOW2, structwindow);
	DDV_MinMaxInt(pDX, structwindow, 0, 999999999);
	DDX_Text(pDX, IDC_ALIGNWINDOW, alignwindow);
	DDV_MinMaxInt(pDX, alignwindow, 0, 9999999);
	DDX_Check(pDX, IDC_SAVECHECK, savefilecheck);
}


BEGIN_MESSAGE_MAP(CDynalign, CFormView)
	ON_MESSAGE(ID_DYNDONE, DynDone)
	//{{AFX_MSG_MAP(CDynalign)
	ON_BN_CLICKED(IDC_SEQUENCEBUTTON, OnSequencebutton)
	ON_BN_CLICKED(IDC_SEQUENCE2BUTTON, OnSequence2button)
	ON_BN_CLICKED(IDC_CTBUTTON, OnCtbutton)
	ON_BN_CLICKED(IDC_CTBUTTON2, OnCtbutton2)
	ON_BN_CLICKED(IDC_ALIGNMENTBUTTON, OnAlignmentbutton)
	ON_BN_CLICKED(IDC_START, OnStart)
	//}}AFX_MSG_MAP
	ON_COMMAND(ID_S1_BP, OnS1Bp)
	ON_COMMAND(ID_S1_CM, OnS1Cm)
	ON_COMMAND(ID_S1_DS, OnS1Ds)
	ON_COMMAND(ID_S1_FMN, OnS1Fmn)
	ON_COMMAND(ID_S1_SS, OnS1Ss)
	ON_COMMAND(ID_S1_PB, OnS1Pb)
	ON_COMMAND(ID_S1_CURRENT, OnS1Current)
	ON_COMMAND(ID_S1_RESET, OnS1Reset)
	ON_COMMAND(ID_S1_SAVE, OnS1Save)
	ON_COMMAND(ID_S1_RESTORE, OnS1Restore)
	ON_COMMAND(ID_S2_BP, OnS2Bp)
	ON_COMMAND(ID_S2_CM, OnS2Cm)
	ON_COMMAND(ID_S2_DS, OnS2Ds)
	ON_COMMAND(ID_S2_FMN, OnS2Fmn)
	ON_COMMAND(ID_S2_SS, OnS2Ss)
	ON_COMMAND(ID_S2_PB, OnS2Pb)
	ON_COMMAND(ID_S2_CURRENT, OnS2Current)
	ON_COMMAND(ID_S2_RESET, OnS2Reset)
	ON_COMMAND(ID_S2_SAVE, OnS2Save)
	ON_COMMAND(ID_S2_RESTORE, OnS2Restore)
	ON_COMMAND(ID_CONSTRAINTSFORALIGNMENT_FORCEALIGNMENT, OnConstraintsforalignmentForcealignment)
	ON_COMMAND(ID_CONSTRAINTSFORALIGNMENT_SHOWCURRENTALIGNMENTCONSTRAINTS, OnConstraintsforalignmentShowcurrentalignmentconstraints)
	ON_COMMAND(ID_CONSTRAINTSFORALIGNMENT_RESETALIGNMENTCONTRAINTS, OnConstraintsforalignmentResetalignmentcontraints)
	ON_COMMAND(ID_CONSTRAINTSFORALIGNMENT_SAVEALIGNMENT, OnConstraintsforalignmentSavealignment)
	ON_COMMAND(ID_CONSTRAINTSFORALIGNMENT_RESTOREALIGNMENT, OnConstraintsforalignmentRestorealignment)
	ON_COMMAND(ID_TEMPERATURE, OnTemperature)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CDynalign diagnostics

#ifdef _DEBUG
void CDynalign::AssertValid() const
{
	CFormView::AssertValid();
}

void CDynalign::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CDynalign message handlers

void CDynalign::OnSequencebutton() 
{
	//Open the standard dialog to pick a sequence
	char *ctname;
	short int i;
	CFileDialog *filedialog;
	
	

	
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Sequence Files (*.seq)|*.seq|Genbank Files|*.gen|Plain Text Files|*.txt||");

	
	filedialog->m_ofn.lpstrInitialDir=GetDynDocument()->startpath;
	if (filedialog->DoModal()==IDOK) {
		
		m_sequencename=(filedialog->GetPathName()).GetBuffer(30);
		i = m_sequencename.GetLength();
		
		ctname = new char[i+4];//allocate enough space so that 
														//three characters can be added 
														//to the name if necessary
		strcpy(ctname,m_sequencename.GetBuffer(10));
		//count the characters to the .
		
		while(i>=0){
			
			if (ctname[i]=='.') break;
			i--;
		}
		if (i==0) i = m_sequencename.GetLength();
		strcpy(ctname+i+1,"ct\0");
		m_ctname=ctname;

		if (m_sequence2name!="") NewSequences();
		
		delete[] ctname;//fix this?
		UpdateData(FALSE);
		
		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(GetDynDocument()->startpath,_MAX_PATH);

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
		strncpy(GetDynDocument()->startpath,path.GetBuffer(1),i);
		*(GetDynDocument()->startpath + i) ='\0';

		
		

	}
	delete filedialog;
	
}

void CDynalign::OnSequence2button() 
{
	//Open the standard dialog to pick a sequence
	char *ctname;
	short int i;
	CFileDialog *filedialog;
	
	

	
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Sequence Files (*.seq)|*.seq|Genbank Files|*.gen|Plain Text Files|*.txt||");

	
	filedialog->m_ofn.lpstrInitialDir=GetDynDocument()->startpath;
	if (filedialog->DoModal()==IDOK) {
		
		m_sequence2name=(filedialog->GetPathName()).GetBuffer(30);
		i = m_sequence2name.GetLength();
		
		ctname = new char[i+4];//allocate enough space so that 
														//three characters can be added 
														//to the name if necessary
		strcpy(ctname,m_sequence2name.GetBuffer(10));
		//count the characters to the .
		
		while(i>=0){
			
			if (ctname[i]=='.') break;
			i--;
		}
		if (i==0) i = m_sequence2name.GetLength();
		strcpy(ctname+i+1,"ct\0");

		m_ctname2=ctname;
		if (m_sequencename!="") {
			NewSequences();



		}

		
		
		delete[] ctname;
		UpdateData(FALSE);
		
		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(GetDynDocument()->startpath,_MAX_PATH);

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
		strncpy(GetDynDocument()->startpath,path.GetBuffer(1),i);
		*(GetDynDocument()->startpath + i) ='\0';

		
		

	}
	delete filedialog;

	
	
}

void CDynalign::OnCtbutton() 
{
	//The user is specifying the CT file name explicitly
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".ct",NULL,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"CT Files (*.ct)|*.ct|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		m_ctname=(filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);		

	}
	delete filedialog;
	
}



void CDynalign::OnCtbutton2() 
{
	//The user is specifying the CT file name explicitly
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".ct",NULL,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"CT Files (*.ct)|*.ct|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		m_ctname2=(filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);		

	}
	delete filedialog;
	
}

void CDynalign::OnAlignmentbutton() 
{
	
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".ct",NULL,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"Alignment Files (*.ali)|*.ali|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		m_ctname2=(filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);		

	}
	delete filedialog;
	
}


CDynDoc *CDynalign::GetDynDocument() {
	
	return ((CDynDoc*) GetDocument());	

}

void CDynalign::NewSequences() {
	//both sequences are defined:
	char *alignmentname,*temp;
	unsigned int i,j;


	alignmentname = new char [m_sequencename.GetLength()+m_sequence2name.GetLength()];
	temp = new char [m_sequence2name.GetLength()+1];

	strcpy(alignmentname,m_sequencename.GetBuffer(0));
	i = strlen(alignmentname);
	while (alignmentname[i]!='.'&&i!=0) i--;

	if (i!=0) alignmentname[i+1]='\0';

	strcpy(temp,m_sequence2name.GetBuffer(0));
	
	j = 0;
	for (i=0;i<=strlen(temp);i++) {
		if (temp[i]=='\\') j = i;

	}

	i = strlen(temp);
	while (temp[i]!='.'&&i!=0) i--;
	
	if (i!=0) temp[i+1] = '\0';
	strcat(alignmentname,temp+j+1);
	strcat(alignmentname,"ali");

	m_alignment = alignmentname;


	GetDynDocument()->ct.openseq(m_sequencename.GetBuffer(0));
	GetDynDocument()->ct2.openseq(m_sequence2name.GetBuffer(0));

	//set M based on the length of sequence 1:
	//this is now disabled because the posterior probs from sequence alignment are better

	//if (GetDynDocument()->ct.numofbases>500) m_M = 15;
	//else if (GetDynDocument()->ct.numofbases>200) m_M = 10;
	//else if (GetDynDocument()->ct.numofbases>75) m_M = 7;
	//else m_M = 4;

	//if (m_M>GetDynDocument()->ct.numofbases) m_M = GetDynDocument()->ct.numofbases;

	if ((GetDynDocument()->ct.GetSequenceLength())>1200) {
   			structwindow=20;
			alignwindow=3;
	}
	else if ((GetDynDocument()->ct.GetSequenceLength())>800) {
   			structwindow=15;
			alignwindow=3; 
	}
	else if ((GetDynDocument()->ct.GetSequenceLength())>500) {
   			structwindow=11;
           alignwindow=3; 
	}
	else if ((GetDynDocument()->ct.GetSequenceLength())>300) {
   			structwindow=7;
			alignwindow=2; 
	}
	else if ((GetDynDocument()->ct.GetSequenceLength())>120) {
   			structwindow=5;
			alignwindow=1;
            
	}
	else if ((GetDynDocument()->ct.GetSequenceLength())>50) {
   			structwindow=3;
			alignwindow=1;
            
	}
	else {
		structwindow=2;
		alignwindow=0;
	}

	maxtracebacks = 20;
	percent = 20;

	delete[] alignmentname;
	delete[] temp;
	UpdateData(FALSE);


}

void CDynalign::OnStart() 
{
	short i,j,k;
	char *temp;
	if (m_sequencename==""||m_sequence2name==""||m_alignment=="") {

		AfxMessageBox( "Please fill all fields by clicking on the button to the left.", 
			MB_OK|MB_ICONHAND);

		return;
	}

	//check to see if the temperature has been changed.
	if (GetDynDocument()->T<310||GetDynDocument()->T>311) {

		//change the temperature from 310.15
		if (GetDynDocument()->newtemp()==0) {
			//if newtemp returned zero, pass a warning to the user
			AfxMessageBox( "An enthalpy data file could not be found!\nTemperature of prediction will revert back to 37 degrees C.\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
				MB_OK|MB_ICONHAND);

		}

	}

	progress = new TProgressDialog();//parent,"Folding the RNA..."
	CRect *rect;
	rect = new CRect(10,40,210,120);
	
	
	progress->Create(NULL,"Dynalign in Progress...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);//(IDD_PROGRESSDIALOG);//
	
	delete rect;
	UpdateData(TRUE);

	dynobject.progress = progress;
	dynobject.ct1 = &GetDynDocument()->ct;
	dynobject.ct2 = &GetDynDocument()->ct2;
	dynobject.data = &GetDynDocument()->data;
    dynobject.parent = this;
	dynobject.ctoutfile = m_ctname.GetBuffer(0);
	dynobject.ctoutfile2 = m_ctname2.GetBuffer(0);
	dynobject.dggap = (short) (m_dggap*10);
	dynobject.aoutfile = m_alignment.GetBuffer(0);
	dynobject.structwindow = structwindow;
	dynobject.alignwindow = alignwindow;
	dynobject.maxtracebacks = maxtracebacks;
	dynobject.percent = percent;
	dynobject.forcealign = forcealign;
	dynobject.datapath = GetDynDocument()->Datapath;
	
	
	dynobject.m = m_M;
	
	if (m_singleinsert)
		dynobject.singleinsert = true;
	else dynobject.singleinsert = false;

	if (savefilecheck) {
		*GetDynDocument()->checksave = TRUE;
		i = m_ctname.GetLength()+m_ctname2.GetLength();
		
		dynobject.savefile = new char[i+4];
		strcpy(dynobject.savefile,m_ctname.GetBuffer(10));
		i = strlen(dynobject.savefile);
		//count the characters to the .
		
		while(i>=0){
			
			if (dynobject.savefile[i]=='.') break;
			i--;
		}
		if (i==0) i = strlen(dynobject.savefile);

		temp = new char[m_ctname2.GetLength()+1];
		strcpy(temp,m_ctname2.GetBuffer(10));
		//start at the right and eat away the path
		for (j=0;(size_t) j<=strlen(temp);j++) {
			if (temp[j]=='\\') k = j+1;	

		}

		strcpy(dynobject.savefile+i+1,temp+k);

		i = strlen(dynobject.savefile);
		//count the characters to the .
		
		while(i>=0){
			
			if (dynobject.savefile[i]=='.') break;
			i--;
		}
		if (i==0) i = strlen(dynobject.savefile);
		strcpy(dynobject.savefile+i+1,"dsv\0");
		delete[] temp;

	}
	else {
		dynobject.savefile=NULL;
		*GetDynDocument()->checksave = FALSE;
	}
	
	AfxBeginThread(DynProc,&dynobject);
	

	
}

LRESULT CDynalign::DynDone(WPARAM wParam, LPARAM lParam) {

	CDialog* finished;
	
	if (savefilecheck) {
		delete[] dynobject.savefile;
	}
	finished = new CDialog(IDD_FINISHED,this);

	if(finished->DoModal()==IDOK) {
		((CRNAstructureApp*) GetDynDocument()->pMainFrame)->Draw(m_ctname.GetBuffer(10));
		((CRNAstructureApp*) GetDynDocument()->pMainFrame)->Draw(m_ctname2.GetBuffer(10));
	}
	

	delete finished;

	progress->SendMessage (WM_CLOSE);
	
	GetDynDocument()->Frame->SendMessage(WM_CLOSE);
	return 0;
}

void CDynalign::OnS1Bp()
{
	CPairDialog *pairdialog;

	pairdialog = new CPairDialog(&(GetDynDocument()->ct),this);

	pairdialog->Create(IDD_PAIR_DIALOG,this);
}

void CDynalign::OnS1Cm()
{
	CSingleDialog *singledialog;

	singledialog = new CSingleDialog(&(GetDynDocument()->ct),this,3);

	singledialog->Create(IDD_SINGLE_DIALOG,this);
}

void CDynalign::OnS1Ds()
{
	CDoubleDialog *doubledialog;

	doubledialog = new CDoubleDialog(&(GetDynDocument()->ct),this);

	doubledialog->Create(IDD_DOUBLE_DIALOG,this);
}

void CDynalign::OnS1Fmn()
{
	CSingleDialog *singledialog;

	singledialog = new CSingleDialog(&(GetDynDocument()->ct),this,2);

	singledialog->Create(IDD_SINGLE_DIALOG,this);
}

void CDynalign::OnS1Ss()
{
	CSingleDialog *singledialog;

	singledialog = new CSingleDialog(&(GetDynDocument()->ct),this);

	singledialog->Create(IDD_SINGLE_DIALOG,this);
}

void CDynalign::OnS1Pb()
{
	CPairDialog *pairdialog;

	pairdialog = new CPairDialog(&(GetDynDocument()->ct),this,true);

	pairdialog->Create(IDD_PAIR_DIALOG,this);
}

void CDynalign::OnS1Current()
{
	currentconstraints(&(GetDynDocument()->ct));
}

void CDynalign::OnS1Reset()
{
	if (AfxMessageBox( "This will erase all constraints.\nContinue?", 
		MB_OKCANCEL   |MB_ICONQUESTION   )==IDOK) {
		//reset all constraints:
		GetDynDocument()->ct.RemoveConstraints();
		//GetDynDocument()->ct.ndbl = 0;
		//GetDynDocument()->ct.npair = 0;
		//GetDynDocument()->ct.nnopair = 0;
		//GetDynDocument()->ct.nmod = 0;
		//GetDynDocument()->ct.ngu=0;
		//GetDynDocument()->ct.nforbid=0;
	}
}

void CDynalign::OnS1Save()
{
	CFileDialog *filedialog;


	filedialog = new CFileDialog(FALSE,".con",NULL,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"Constraint Files (*.con)|*.con|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {


		outputconstraints(filedialog->GetPathName().GetBuffer(0),&(GetDynDocument()->ct));

	

	}
	delete filedialog;
}

void CDynalign::OnS1Restore()
{
	CFileDialog *filedialog;
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Constraint Files (*.con)|*.con|All Files|*.*||");

	
	
	if (filedialog->DoModal()==IDOK) {
		readconstraints(filedialog->GetPathName().GetBuffer(0),&(GetDynDocument()->ct));
		

	}
	delete filedialog;
}

void CDynalign::OnS2Bp()
{
	CPairDialog *pairdialog;

	pairdialog = new CPairDialog(&(GetDynDocument()->ct2),this);

	pairdialog->Create(IDD_PAIR_DIALOG,this);
}

void CDynalign::OnS2Cm()
{
	CSingleDialog *singledialog;

	singledialog = new CSingleDialog(&(GetDynDocument()->ct2),this,3);

	singledialog->Create(IDD_SINGLE_DIALOG,this);
}

void CDynalign::OnS2Ds()
{
	CDoubleDialog *doubledialog;

	doubledialog = new CDoubleDialog(&(GetDynDocument()->ct2),this);

	doubledialog->Create(IDD_DOUBLE_DIALOG,this);
}

void CDynalign::OnS2Fmn()
{
	CSingleDialog *singledialog;

	singledialog = new CSingleDialog(&(GetDynDocument()->ct2),this,2);

	singledialog->Create(IDD_SINGLE_DIALOG,this);
}

void CDynalign::OnS2Ss()
{
	CSingleDialog *singledialog;

	singledialog = new CSingleDialog(&(GetDynDocument()->ct2),this);

	singledialog->Create(IDD_SINGLE_DIALOG,this);
}

void CDynalign::OnS2Pb()
{
	CPairDialog *pairdialog;

	pairdialog = new CPairDialog(&(GetDynDocument()->ct2),this,true);

	pairdialog->Create(IDD_PAIR_DIALOG,this);
}

void CDynalign::OnS2Current()
{
	currentconstraints(&(GetDynDocument()->ct2));
}

void CDynalign::OnS2Reset()
{
	if (AfxMessageBox( "This will erase all constraints.\nContinue?", 
		MB_OKCANCEL   |MB_ICONQUESTION   )==IDOK) {
		//reset all constraints:
			GetDynDocument()->ct2.RemoveConstraints();
		//GetDynDocument()->ct2.ndbl = 0;
		//GetDynDocument()->ct2.npair = 0;
		//GetDynDocument()->ct2.nnopair = 0;
		//GetDynDocument()->ct2.nmod = 0;
		//GetDynDocument()->ct2.ngu=0;
		//GetDynDocument()->ct2.nforbid=0;
	}
}

void CDynalign::OnS2Save()
{
	CFileDialog *filedialog;


	filedialog = new CFileDialog(FALSE,".con",NULL,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"Constraint Files (*.con)|*.con|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {


		outputconstraints(filedialog->GetPathName().GetBuffer(0),&(GetDynDocument()->ct2));

	

	}
	delete filedialog;
}

void CDynalign::OnS2Restore()
{
	CFileDialog *filedialog;
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Constraint Files (*.con)|*.con|All Files|*.*||");

	
	
	if (filedialog->DoModal()==IDOK) {
		readconstraints(filedialog->GetPathName().GetBuffer(0),&(GetDynDocument()->ct2));
		

	}
	delete filedialog;
}

void CDynalign::OnConstraintsforalignmentForcealignment()
{
	CAlignForceDialog *alignforcedialog;
	int index;
	

	if (forcealign==NULL) {
		//allocate forcealign
		forcealign=new short *[2];
		
		forcealign[0]=new short [GetDynDocument()->ct.GetSequenceLength()+1];
		forcealign[1]=new short [GetDynDocument()->ct2.GetSequenceLength()+1];
		for (index=1;index<=GetDynDocument()->ct.GetSequenceLength();index++) forcealign[0][index]=0;
		for (index=1;index<=GetDynDocument()->ct2.GetSequenceLength();index++) forcealign[1][index]=0;
		


	}

	alignforcedialog = new CAlignForceDialog(forcealign,GetDynDocument()->ct.GetSequenceLength(),GetDynDocument()->ct2.GetSequenceLength(),this);

	alignforcedialog->Create(IDD_ALIGNDIALOG,this);
}

void CDynalign::OnConstraintsforalignmentShowcurrentalignmentconstraints()
{
	CString message;
	int index;
	char digit[100];

	if (forcealign==NULL) {

		message="There are no alignment constraints.";
	}
	else {
		message="";
		for (index=1;index<=GetDynDocument()->ct.GetSequenceLength();index++) {

			if (forcealign[0][index]>0) {
				message+="Nucleotide #";
				sprintf(digit,"%i",index);
				message+=digit;
				message+=" in sequence 1 aligned to #";
				sprintf(digit,"%i",forcealign[0][index]);
				message+=digit;
				message+=" in sequence 2.\n\n";


			}

		}

	}


	AfxMessageBox( message, 
			MB_OK|MB_ICONINFORMATION   );
}

void CDynalign::OnConstraintsforalignmentResetalignmentcontraints()
{
	if (AfxMessageBox( "This will erase all alignmnet constraints.\nContinue?", 
		MB_OKCANCEL   |MB_ICONQUESTION   )==IDOK) {
		//reset all constraints:
		if (forcealign!=NULL) {
			delete[] forcealign[0];
			delete[] forcealign[1];
			delete[] forcealign;
			forcealign=NULL;

		}
	}
}

void CDynalign::OnConstraintsforalignmentSavealignment()
{
	CFileDialog *filedialog;
	ofstream out;
	int index;


	filedialog = new CFileDialog(FALSE,".con",NULL,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"Alignment Constraint Files (*.acon)|*.acon|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {


		out.open(filedialog->GetPathName().GetBuffer(0));
		out << GetDynDocument()->ct.GetSequenceLength() << "\n";
		out << GetDynDocument()->ct2.GetSequenceLength() << "\n";
		for (index=1;index<=GetDynDocument()->ct.GetSequenceLength();index++) {

			out << forcealign[0][index]<<"\n";
		}
		for (index=1;index<=GetDynDocument()->ct2.GetSequenceLength();index++) {

			out << forcealign[1][index]<<"\n";
		}

		out.close();

	}
	delete filedialog;
}

void CDynalign::OnConstraintsforalignmentRestorealignment()
{
	CFileDialog *filedialog;
	ifstream in;
	int temp;
	int index;

	if (AfxMessageBox( "This will erase all current alignmnet constraints.\nContinue?", 
		MB_OKCANCEL   |MB_ICONQUESTION   )==IDOK) {
		//reset all constraints:
		if (forcealign!=NULL) {
			delete[] forcealign[0];
			delete[] forcealign[1];
			delete[] forcealign;
			forcealign=NULL;

		}
		

		
	
		filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
			"Constraint Files (*.acon)|*.acon|All Files|*.*||");

	
	
		if (filedialog->DoModal()==IDOK) {
			forcealign=new short *[2];
		
			forcealign[0]=new short [GetDynDocument()->ct.GetSequenceLength()+1];
			forcealign[1]=new short [GetDynDocument()->ct2.GetSequenceLength()+1];
			for (index=1;index<=GetDynDocument()->ct.GetSequenceLength();index++) forcealign[0][index]=0;
			for (index=1;index<=GetDynDocument()->ct2.GetSequenceLength();index++) forcealign[1][index]=0;

			in.open(filedialog->GetPathName().GetBuffer(0));
			in >> temp;
			if (temp!=GetDynDocument()->ct.GetSequenceLength()) {
				AfxMessageBox( "The alignment save file was created for a different length sequence 1.\nNo constraints will be loaded.", 
					MB_OK|MB_ICONINFORMATION   );
				in.close();
				delete filedialog;
				return;

			}
			in >> temp;
			if (temp!=GetDynDocument()->ct2.GetSequenceLength()) {
				AfxMessageBox( "The alignment save file was created for a different length sequence 2.\nNo constraints will be loaded.", 
					MB_OK|MB_ICONINFORMATION   );
				in.close();
				
				delete filedialog;
				return;

			}

			
			for (index=1;index<=GetDynDocument()->ct.GetSequenceLength();index++) {

				in >> forcealign[0][index];
			}
			for (index=1;index<=GetDynDocument()->ct2.GetSequenceLength();index++) {

				in >> forcealign[1][index];
			}

			in.close();
		

		}
		delete filedialog;
		

	}
}

void CDynalign::OnTemperature()
{
	//Allow the user to specify a new temperature.
	CTemp_Dialog *temp;
	temp=new CTemp_Dialog(&(GetDynDocument()->T));

	temp->DoModal();
	delete temp;
	


}


