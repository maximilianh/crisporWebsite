// BiFoldView.cpp : implementation file
//
#include <string>
#include "stdafx.h"
#include "RNAstructure.h"
#include "BiFoldView.h"
#include "advanced.h"
#include "temp_Dialog.h"

#include <direct.h>



#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CBiFoldView dialog
IMPLEMENT_DYNCREATE(CBiFoldView, CFormView)

CBiFoldView::CBiFoldView(CWnd* pParent /*=NULL*/)
	: CFormView(CBiFoldView::IDD)
	 
{
	//{{AFX_DATA_INIT(CBiFoldView)
	m_ctname = _T("");
	m_sequence1 = _T("");
	m_sequence2 = _T("");
	m_window = 0;
	m_percent = 50;
	m_number = 20;
	m_save = FALSE;
	//}}AFX_DATA_INIT

	dynwindow= 0;
	dynnumber=750;
	dynpercent=20;

	forbidunimolecular = false;

	started = false;//used to track whether the calculation has started
}

//Protected constructor
//CBiFoldView::CBiFoldView()
//	: CFormView(CBiFoldView::IDD)
//{
	
//}

void CBiFoldView::OnInitialUpdate() {

	if (GetBiFoldDocument()->savefile!=NULL) {
		if (*GetBiFoldDocument()->savefile) {
			m_save = TRUE;
			UpdateData(FALSE);

		}
	}

	ResizeParentToFit(FALSE);
	ResizeParentToFit();
}

void CBiFoldView::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CBiFoldView)
	DDX_Text(pDX, IDC_CTNAME, m_ctname);
	DDX_Text(pDX, IDC_SEQUENCENAME, m_sequence1);
	DDX_Text(pDX, IDC_SEQUENCE2NAME, m_sequence2);
	DDX_Text(pDX, IDC_WINDOW, m_window);
	DDV_MinMaxInt(pDX, m_window, 0, 16000);
	DDX_Text(pDX, IDC_PERCENT, m_percent);
	DDV_MinMaxInt(pDX, m_percent, 0, 99);
	DDX_Text(pDX, IDC_NUMBER, m_number);
	//}}AFX_DATA_MAP
	DDX_Check(pDX, IDC_SAVECHECK, m_save);
}



/////////////////////////////////////////////////////////////////////////////
// CBiFoldView message handlers

void CBiFoldView::OnUpdate(CView*, LPARAM, CObject*)
{

	ResizeParentToFit(FALSE);
	ResizeParentToFit();
	

}

void CBiFoldView::OnSequencebutton() 
{
	
	
	CFileDialog *filedialog;
	
	

	
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Sequence Files (*.seq)|*.seq||");

	
	filedialog->m_ofn.lpstrInitialDir=GetBiFoldDocument()->startpath;
	if (filedialog->DoModal()==IDOK) {
		//strcpy(m_sequencename.GetBuffer(10),(filedialog->GetPathName()).GetBuffer(0));
		m_sequence1=(filedialog->GetPathName()).GetBuffer(30);
		
		
		
		UpdateData(FALSE);
		
		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(GetBiFoldDocument()->startpath,_MAX_PATH);
		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(GetFoldDocument()->startpath,_MAX_PATH);
		CString path;
		path = filedialog->GetPathName();
		int i = path.GetLength();
		while(i>=0){
			
			if (path[i]=='\\') break;
			i--;
		}
		if (i>_MAX_PATH) i = _MAX_PATH;
		strncpy(GetBiFoldDocument()->startpath,path.GetBuffer(1),i);
		*(GetBiFoldDocument()->startpath + i) ='\0';
		
		

	}
	delete filedialog;

	
}


void CBiFoldView::OnSequence2button() 
{
	CFileDialog *filedialog;
	char *seqname1,*seqname2;
	int i;

	
	

	
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Sequence Files (*.seq)|*.seq||");

	
	filedialog->m_ofn.lpstrInitialDir=GetBiFoldDocument()->startpath;
	if (filedialog->DoModal()==IDOK) {
		//strcpy(m_sequencename.GetBuffer(10),(filedialog->GetPathName()).GetBuffer(0));
		m_sequence2=(filedialog->GetPathName()).GetBuffer(30);
		
		
		
		
		
		
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
		strncpy(GetBiFoldDocument()->startpath,path.GetBuffer(1),i);
		*(GetBiFoldDocument()->startpath + i) ='\0';



		if (m_sequence1!="") {
			i = m_sequence1.GetLength();
			seqname1= new char [(i)+1];
			i = m_sequence2.GetLength();
			seqname2= new char [(i)+1];

			strcpy(seqname1,m_sequence1.GetBuffer(1));
			strcpy(seqname2,m_sequence2.GetBuffer(2));

			i = strlen(seqname1);
			while(i>=0){
			
				if (seqname1[i]=='.') break;
				i--;
			}
			if (i!=0) seqname1[i] = '\0';
			
			
			
			
			
			i = strlen(seqname2);
			while(i>=0){
			
				if (seqname2[i]=='.') break;
				i--;
			}
			if (i!=0) seqname2[i] = '\0';

			//get rid of the path from seqname2
			i = strlen(seqname2);
			while (i>=0) {
				if (seqname2[i]=='\\') break;
				i--;

			}
			strcpy(seqname2,seqname2+i+1);
			
			m_ctname+=seqname1;
			m_ctname+="_";
			m_ctname+=seqname2;
			m_ctname+=".ct";

			delete[] seqname1;
			delete[] seqname2;
		}
		
		UpdateData(FALSE);

	}
	delete filedialog;

	
}
void CBiFoldView::OnCtbutton() 
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

CBiFoldDoc *CBiFoldView::GetBiFoldDocument() {
	
	return ((CBiFoldDoc*) GetDocument());	

}

void CBiFoldView::OnStart() 
{
	int i,j;

	//a trap to make sure the calculation is started more than once
	if (started) return;
	started=true;

	UpdateData(TRUE);

	
	if(m_sequence1=="") {
		MessageBox("Please specify a sequence name for the first sequence.");
		return;
	}
	else if(m_sequence2=="") {
		MessageBox("Please specify a sequence name for the second sequence./nUse unimolecular folding for a single sequence.");
		return;
	}

	

	if(GetBiFoldDocument()->ct.openseq(m_sequence1.GetBuffer(10))==0||
		(ct2.openseq (m_sequence2.GetBuffer(10)))==0) {

		AfxMessageBox( "One of the sequences is unrecognizable.", 
			MB_OK|MB_ICONSTOP);

		return;

	}

    //if (dynnumber<m_number) dynnumber = m_number;
    //if (dynpercent<m_percent) dynpercent = m_percent;
    //if (dynwindow>m_window) dynwindow = m_window;

	string stringname;
	stringname = GetBiFoldDocument()->ct.GetSequenceLabel();
	//strcpy(ct3.ctlabel[1],GetBiFoldDocument()->ct.ctlabel[1]);
	//remove the new line at the end of ct.ctlabel[1]
	if (stringname[stringname.size()-1]=='\n') {
		stringname.erase(stringname.size()-1,1);		
	}

	stringname+="_";
	stringname+=ct2.GetSequenceLabel();

	//i = strlen(ct3.ctlabel[1]);
	//ct3.ctlabel[1][i-1]='\0';

	//strcat(ct3.ctlabel[1],"_");
	//strcat(ct3.ctlabel[1],ct2.ctlabel[1]);

	ct3.SetSequenceLabel(stringname);

	//prepare ct3 with both sequences:
	/*ct3.numofbases=GetBiFoldDocument()->ct.numofbases+ct2.numofbases+3;
	ct3.allocate(GetBiFoldDocument()->ct.numofbases+ct2.numofbases+3);
	for (i=1;i<=GetBiFoldDocument()->ct.numofbases;i++) {
		ct3.numseq[i] = GetBiFoldDocument()->ct.numseq[i];
		ct3.nucs[i] = GetBiFoldDocument()->ct.nucs[i];
		ct3.hnumber[i] = GetBiFoldDocument()->ct.hnumber[i];

	}
	
	for (i=1;i<=ct2.numofbases;i++) {
		ct3.numseq[i+GetBiFoldDocument()->ct.numofbases+3] = ct2.numseq[i];
		ct3.nucs[i+GetBiFoldDocument()->ct.numofbases+3] = ct2.nucs[i];
		ct3.hnumber[i+GetBiFoldDocument()->ct.numofbases+3] = ct2.hnumber[i];

	} 	
      
   
   ct3.numseq[GetBiFoldDocument()->ct.numofbases+1] = 5;
   ct3.numseq[GetBiFoldDocument()->ct.numofbases+2] = 5;
   ct3.numseq[GetBiFoldDocument()->ct.numofbases+3] = 5;

   ct3.nucs[GetBiFoldDocument()->ct.numofbases+1] = 'I';
   ct3.nucs[GetBiFoldDocument()->ct.numofbases+2] = 'I';
   ct3.nucs[GetBiFoldDocument()->ct.numofbases+3] = 'I';

   ct3.hnumber[GetBiFoldDocument()->ct.numofbases+1] = 0;
   ct3.hnumber[GetBiFoldDocument()->ct.numofbases+2] = 0;
   ct3.hnumber[GetBiFoldDocument()->ct.numofbases+3] = 0;


   ct3.inter[0] = GetBiFoldDocument()->ct.numofbases+1;
   ct3.inter[1] = GetBiFoldDocument()->ct.numofbases+2;
   ct3.inter[2] = GetBiFoldDocument()->ct.numofbases+3;

   ct3.intermolecular = true;*/

	ct3.allocate(GetBiFoldDocument()->ct.GetSequenceLength()+ct2.GetSequenceLength()+3);

	for (i=1;i<=GetBiFoldDocument()->ct.GetSequenceLength();i++) {
		ct3.numseq[i] = GetBiFoldDocument()->ct.numseq[i];
		ct3.nucs[i] = GetBiFoldDocument()->ct.nucs[i];
		ct3.hnumber[i] = GetBiFoldDocument()->ct.hnumber[i];

	}
	
	for (i=1;i<=ct2.GetSequenceLength();i++) {
		ct3.numseq[i+GetBiFoldDocument()->ct.GetSequenceLength()+3] = ct2.numseq[i];
		ct3.nucs[i+GetBiFoldDocument()->ct.GetSequenceLength()+3] = ct2.nucs[i];
		ct3.hnumber[i+GetBiFoldDocument()->ct.GetSequenceLength()+3] = ct2.hnumber[i];

	} 	
      
   
   ct3.numseq[GetBiFoldDocument()->ct.GetSequenceLength()+1] = 5;
	ct3.numseq[GetBiFoldDocument()->ct.GetSequenceLength()+2] = 5;
   ct3.numseq[GetBiFoldDocument()->ct.GetSequenceLength()+3] = 5;

   ct3.nucs[GetBiFoldDocument()->ct.GetSequenceLength()+1] = 'I';
   ct3.nucs[GetBiFoldDocument()->ct.GetSequenceLength()+2] = 'I';
   ct3.nucs[GetBiFoldDocument()->ct.GetSequenceLength()+3] = 'I';

   ct3.hnumber[GetBiFoldDocument()->ct.GetSequenceLength()+1] = 0;
   ct3.hnumber[GetBiFoldDocument()->ct.GetSequenceLength()+2] = 0;
   ct3.hnumber[GetBiFoldDocument()->ct.GetSequenceLength()+3] = 0;


   ct3.inter[0] = GetBiFoldDocument()->ct.GetSequenceLength()+1;
   ct3.inter[1] = GetBiFoldDocument()->ct.GetSequenceLength()+2;
   ct3.inter[2] = GetBiFoldDocument()->ct.GetSequenceLength()+3;

   ct3.intermolecular = true;


   //Also copy information about nucleotides that must be single stranded
		//(These were entered as lowercase by the user.)

   for (i=0;i<GetBiFoldDocument()->ct.GetNumberofSingles();i++) {
		ct3.AddSingle(GetBiFoldDocument()->ct.GetSingle(i));

   }
   for (i=0;i<ct2.GetNumberofSingles();i++) {
	   ct3.AddSingle(ct2.GetSingle(i)+GetBiFoldDocument()->ct.GetSequenceLength()+3);
		//ct3.nnopair++;
		//ct3.nopair[ct3.nnopair] = ct2.nopair[i]+GetBiFoldDocument()->ct.numofbases+3;

   }


   if (forbidunimolecular) {
		//forbid unimolecular pairs
	   ct3.allocatetem();
	   for (i=1;i<GetBiFoldDocument()->ct.GetSequenceLength();i++) {
		   for (j=i+1;j<=GetBiFoldDocument()->ct.GetSequenceLength();j++) {
				ct3.tem[j][i]=false;
		   }

	   }
	   for (i=GetBiFoldDocument()->ct.GetSequenceLength()+3;i<ct3.GetSequenceLength();i++) {
		   for (j=i+1;j<=ct3.GetSequenceLength();j++) {
				ct3.tem[j][i]=false;
		   }

	   }
   }

   //check to see if the temperature has been changed.
	if (GetBiFoldDocument()->T<310||GetBiFoldDocument()->T>311) {

		//change the temperature from 310.15
		if (GetBiFoldDocument()->newtemp()==0) {
			//if newtemp returned zero, pass a warning to the user
			AfxMessageBox( "An enthalpy data file could not be found!\nTemperature of prediction will revert back to 37 degrees C.\nData files can be downloaded on the Mathews Lab website,\nhttp://rna.urmc.rochester.edu.", 
				MB_OK|MB_ICONHAND);

		}

	}


	//TProgressDialog *progress;
	progress = new TProgressDialog();//parent,"Folding the RNA..."
	CRect *rect;
	rect = new CRect(10,40,210,120);
	
	if (GetBiFoldDocument()->ISRNA)
	progress->Create(NULL,"Folding the RNA...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);
	else progress->Create(NULL,"Folding the DNA...",WS_CHILD|WS_VISIBLE|WS_CAPTION,*rect);
	delete rect;


	
	
	foldobject.parent=this;
	foldobject.ct=&ct3;
	foldobject.data=&GetBiFoldDocument()->data;
	foldobject.dynpercent = m_percent;
	foldobject.dynwindow = m_window;
	//foldobject.m_number = m_number;
	//foldobject.m_percent = m_percent;
	//foldobject.m_window = m_window;
	foldobject.progress = progress;
	foldobject.ctoutfile = m_ctname.GetBuffer(10);
	foldobject.dynnumber = m_number;
	foldobject.subfold=false;//No "subfolding" allowed for bimolecular. 
							//This is a feature for single strands only. 


	//The maximum internal loop size is fixed at 30 for now.
	foldobject.maximuminternalloopsize=30;

	if (m_save) {
		i = m_ctname.GetLength();
		
		foldobject.savefile = new char[i+4];
		strcpy(foldobject.savefile,m_ctname.GetBuffer(10));
		//count the characters to the .
		
		while(i>=0){
			
			if (foldobject.savefile[i]=='.') break;
			i--;
		}
		if (i==0) i = strlen(foldobject.savefile);
		strcpy(foldobject.savefile+i+1,"sav\0");
		if (GetBiFoldDocument()->savefile !=NULL)
			*GetBiFoldDocument()->savefile = true; 

	}
	else {
		foldobject.savefile = 0;
		if (GetBiFoldDocument()->savefile !=NULL)
			*GetBiFoldDocument()->savefile = false; 
	}
	


	



	
	AfxBeginThread(FoldProc,&foldobject);


	

}


LRESULT CBiFoldView::DoneFolding(WPARAM wParam, LPARAM lParam) {
	CDialog* finished;
	

	finished = new CDialog(IDD_FINISHED,this);

	if(finished->DoModal()==IDOK) ((CRNAstructureApp*) GetBiFoldDocument()->pMainFrame)->Draw(m_ctname.GetBuffer(10));
	

	delete finished;

	


	progress->SendMessage (WM_CLOSE);
	
	GetBiFoldDocument()->Frame->SendMessage(WM_CLOSE);
	
	return 0;	


}

/*CFoldDoc *CBiFoldView::GetFoldDocument() {
	
	return ((CFoldDoc*) GetDocument());	

}*/

BEGIN_MESSAGE_MAP(CBiFoldView, CFormView)
	
	//{{AFX_MSG_MAP(CBiFoldView)
	ON_BN_CLICKED(IDC_SEQUENCEBUTTON, OnSequencebutton)
	ON_BN_CLICKED(IDC_CTBUTTON, OnCtbutton)
	ON_BN_CLICKED(IDC_SEQUENCE2BUTTON, OnSequence2button)
	ON_BN_CLICKED(IDC_START, OnStart)
	//}}AFX_MSG_MAP
	ON_MESSAGE(ID_FOLDDONE, DoneFolding)
	ON_COMMAND(ID_FORCE_FORBIDUNIMOLECULARPAIRS, OnForceForbidunimolecularpairs)
	ON_WM_SETFOCUS()
	ON_COMMAND(ID_TEMPERATURE, OnTemperature)
END_MESSAGE_MAP()




void CBiFoldView::OnForceForbidunimolecularpairs()
{
	//User has chosen the forbid unimolecular pair menu option
	CMenu* menu = GetBiFoldDocument()->menuframe->GetMenu( );	

	forbidunimolecular = !forbidunimolecular;

	if (forbidunimolecular) menu->CheckMenuItem(ID_FORCE_FORBIDUNIMOLECULARPAIRS,MF_CHECKED);
	else menu->CheckMenuItem(ID_FORCE_FORBIDUNIMOLECULARPAIRS,MF_UNCHECKED); 

	
}

void CBiFoldView::OnSetFocus(CWnd* pOldWnd)
{
	CFormView::OnSetFocus(pOldWnd);

	//In case the user has two bifoldview windows open, the menu displayed must show the correct checkmark
	//state -- this code guarantees that
	CMenu* menu = GetBiFoldDocument()->menuframe->GetMenu( );
	if (forbidunimolecular) menu->CheckMenuItem(ID_FORCE_FORBIDUNIMOLECULARPAIRS,MF_CHECKED);
	else menu->CheckMenuItem(ID_FORCE_FORBIDUNIMOLECULARPAIRS,MF_UNCHECKED);

	
}

void CBiFoldView::OnTemperature()
{
	//Allow the user to specify a new temperature.
	CTemp_Dialog *temp;
	temp=new CTemp_Dialog(&(GetBiFoldDocument()->T));

	temp->DoModal();
	delete temp;
	


}
