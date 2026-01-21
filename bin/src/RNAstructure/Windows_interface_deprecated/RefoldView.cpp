// RefoldView.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "RefoldView.h"
#include "ReFoldDoc.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#include <direct.h>
#include "../src/outputconstraints.h"

#include <iostream>
using namespace std;

/////////////////////////////////////////////////////////////////////////////
// CRefoldView

IMPLEMENT_DYNCREATE(CRefoldView, CFormView)

CRefoldView::CRefoldView()
	: CFormView(CRefoldView::IDD)
{
	//{{AFX_DATA_INIT(CRefoldView)
	m_ctname = _T("");
	m_number = 20;
	m_percent = 10;
	m_savename = _T("");
	m_window = 0;
	//}}AFX_DATA_INIT
}



CRefoldView::~CRefoldView()
{
}

void CRefoldView::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CRefoldView)
	DDX_Text(pDX, IDC_CTNAME, m_ctname);
	DDX_Text(pDX, IDC_NUMBER, m_number);
	DDX_Text(pDX, IDC_PERCENT, m_percent);
	DDV_MinMaxInt(pDX, m_percent, 0, 100);
	DDX_Text(pDX, IDC_SEQUENCENAME, m_savename);
	DDX_Text(pDX, IDC_WINDOW, m_window);
	DDV_MinMaxInt(pDX, m_window, 0, 16000);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CRefoldView, CFormView)
	//{{AFX_MSG_MAP(CRefoldView)
	ON_BN_CLICKED(IDC_SAVEBUTTON, OnSavebutton)
	ON_BN_CLICKED(IDC_CTBUTTON, OnCtbutton)
	ON_BN_CLICKED(IDC_START, OnStart)
	ON_COMMAND(ID_FORCE_CURRENT, OnForceCurrent)
	ON_COMMAND(ID_FORCE_SAVECONSTRAINTS, OnForceSaveconstraints)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CRefoldView diagnostics

#ifdef _DEBUG
void CRefoldView::AssertValid() const
{
	CFormView::AssertValid();
}

void CRefoldView::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CRefoldView message handlers

void CRefoldView::OnInitialUpdate() 
{
	CFormView::OnInitialUpdate();
	
	ResizeParentToFit();
	
}

void CRefoldView::OnSavebutton() 
{
		//Open the standard dialog to pick a sequence
	char *ctname;
	short int i;//needs to be short int length for reading sequence length
	CFileDialog *filedialog;
	short vers;

		
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Save Files (*.sav)|*.sav||");

	
	filedialog->m_ofn.lpstrInitialDir=((CReFoldDoc*) GetDocument())->startingpath;
	if (filedialog->DoModal()==IDOK) {
		//strcpy(m_sequencename.GetBuffer(10),(filedialog->GetPathName()).GetBuffer(0));

		//get the standard window sizes to fill in
		ifstream sav((filedialog->GetPathName()).GetBuffer(30),ios::binary);

		sav.read((char* ) &vers,sizeof(vers));//read the version of the save file
		
		if (vers!=safiversion) {
			//Wrong version!
			 AfxMessageBox( "This savefile was created with a different version of RNAstructure and cannot be used.", 
			MB_OK|MB_ICONINFORMATION   );
			
		}

		else {


			sav.read((char* ) &i,sizeof(i));
			

			if (i>1200) {
   				m_window=20;
			}
			else if (i>800) {
   				m_window=15;
	            
			}
			else if (i>500) {
   				m_window=11;
	            
			}
			else if (i>300) {
   				m_window=7;
	            
			}
			else if (i>120) {
   				m_window=5;
	            
			}
			else if (i>50) {
   				m_window=3;
	            
			}
			else m_window=2;
			
			m_savename=(filedialog->GetPathName()).GetBuffer(30);
			i = m_savename.GetLength();
			
			ctname = new char[i+4];//allocate enough space so that 
															//three characters can be added 
															//to the name if necessary
			strcpy(ctname,m_savename.GetBuffer(10));
			//count the characters to the .
			
			while(i>=0){
				
				if (ctname[i]=='.') break;
				i--;
			}
			if (i==0) i = m_savename.GetLength();
			strcpy(ctname+i+1,"ct\0");
			m_ctname=ctname;
			
			delete[] ctname;//fix this?
			
			
			//now store the path in Startpath so that the program can start here next time:
			_getcwd(((CReFoldDoc*) GetDocument())->startingpath,_MAX_PATH);
			UpdateData(FALSE);

		}
		sav.close();	
	}
	
	delete filedialog;

	
	
}

void CRefoldView::OnCtbutton() 
{
		//The user is specifying the CT file name explicitly
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".ct",m_ctname,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"CT Files (*.ct)|*.ct|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		m_ctname=(filedialog->GetPathName()).GetBuffer(0);
		UpdateData(FALSE);		

	}
	delete filedialog;
	
}

void CRefoldView::OnStart() 
{
	CDialog* finished;
	UpdateData(TRUE);
	short vers;

	if (m_savename=="") {

		MessageBox("Please specify a save file name.");
		return;

	}

	ifstream sav(m_savename.GetBuffer(30),ios::binary);
	read(&sav,&vers);
	if (vers!=safiversion) {
		//This means that the save file was created with a different code version

		//don't read the file for safety's sake
		MessageBox("This save file is from a different RNAstructure version.");
		return;

	}
	sav.close();
	
	opensav(m_savename.GetBuffer(30),&(((CReFoldDoc*) GetDocument())->ct),m_number,m_percent,m_window);
	(((CReFoldDoc*) GetDocument())->ct).ctout(m_ctname.GetBuffer(30));
	finished = new CDialog(IDD_FINISHED,this);

	if(finished->DoModal()==IDOK) ((CRNAstructureApp*) ((CReFoldDoc*) GetDocument())->pMainFrame)->Draw(m_ctname.GetBuffer(10));
	

	delete finished;
	
	((CReFoldDoc*) GetDocument())->Frame->SendMessage(WM_CLOSE);
}

void CRefoldView::OnForceCurrent() 
{
	char message[5000],temp[5];
	int it,jt,i;

	strcpy (message,"");

	if (m_savename!="") {
		//open the file to get the folding constraints, but do no tracebacks:

		opensav(m_savename.GetBuffer(30),&(((CReFoldDoc*) GetDocument())->ct),0,0,0);


	}



	if (((CReFoldDoc*) GetDocument())->ct.GetNumberofDoubles()>0) {
   	strcat (message,"Bases forced double:\n");
      for (it=0;it<((CReFoldDoc*) GetDocument())->ct.GetNumberofDoubles();it++) {
      	itoa(((CReFoldDoc*) GetDocument())->ct.GetDouble(it),temp,10);
         strcat(message,"   ");
         strcat(message,temp);

         if (it%5==0) strcat(message,"\n");
		 else strcat(message,", ");
      }
	  strcat(message,"\n");
   }


   if (((CReFoldDoc*) GetDocument())->ct.GetNumberofPairs()>0) {
   	strcat(message,"Forced Base Pairs:\n");
   	for (it = 0; it<((CReFoldDoc*) GetDocument())->ct.GetNumberofPairs();it++) {

 		itoa(((CReFoldDoc*) GetDocument())->ct.GetPair5(it),temp,10);
         strcat(message,"  ");
         strcat(message,temp);
         strcat(message," - ");
         itoa(((CReFoldDoc*) GetDocument())->ct.GetPair3(it),temp,10);
         strcat(message,temp);
         if (it%5==0) strcat(message,"\n");
		 else strcat(message,", ");
      }
		strcat(message,"\n");
   }

   if (((CReFoldDoc*) GetDocument())->ct.GetNumberofSingles()>0) {
   	strcat (message,"Bases forced single:\n");
      for (it=0;it<((CReFoldDoc*) GetDocument())->ct.GetNumberofSingles();it++) {
      	itoa(((CReFoldDoc*) GetDocument())->ct.GetSingle(it),temp,10);
         strcat(message,"   ");
         strcat(message,temp);
         if (it%5==0) strcat(message,"\n");
		 else strcat(message,", ");
      }
	  strcat(message,"\n");
   }

   if (((CReFoldDoc*) GetDocument())->ct.GetNumberofGU()>0) {
   	strcat (message,"U's in GU pairs:\n");
      for (it=0;it<((CReFoldDoc*) GetDocument())->ct.GetNumberofGU();it++) {
      	itoa(((CReFoldDoc*) GetDocument())->ct.GetGUpair(it),temp,10);
         strcat(message,"   ");
         strcat(message,temp);
         if (((it+1)%5==0)) strcat(message,"\n");
		 else strcat(message,", ");
      }
	  strcat(message,"\n");
   }

   if (((CReFoldDoc*) GetDocument())->ct.GetNumberofModified()>0) {
   	strcat (message,"Bases accessible to chemical modification:\n");
      for (it=0;it<((CReFoldDoc*) GetDocument())->ct.GetNumberofModified();it++) {
      	itoa(((CReFoldDoc*) GetDocument())->ct.GetModified(it),temp,10);
         strcat(message,"   ");
         strcat(message,temp);
         if (it%5==0) strcat(message,"\n");
		 else strcat(message,", ");
      }
	  strcat(message,"\n");
   }

   

   if (((CReFoldDoc*) GetDocument())->ct.neighbors[0][0]>0) {
	   strcat (message,"Paired neighbors:");
      
      for (it=0;it<100&&((CReFoldDoc*) GetDocument())->ct.neighbors[it][0]>0;it++) {
		  for (jt=0;jt<25&&((CReFoldDoc*) GetDocument())->ct.neighbors[it][jt]>0;jt++) {
			strcat(message,tobase(((CReFoldDoc*) GetDocument())->ct.neighbors[it][jt]));
		  }
		  strcat(message,"\n");
      }
	  strcat(message,"\n");
		 
   }

   if (((CReFoldDoc*) GetDocument())->ct.min_gu>0) {
	  strcat (message,"Minimum # of G-U pairs=");
      
      	itoa(((CReFoldDoc*) GetDocument())->ct.min_gu,temp,10);
         strcat(message," ");
         strcat(message,temp);
         strcat(message,"\n");
		 
   }

   if (((CReFoldDoc*) GetDocument())->ct.min_g_or_u>0) {
	  strcat (message,"Minimum # of G's and U's in pairs=");
      
      	itoa(((CReFoldDoc*) GetDocument())->ct.min_g_or_u,temp,10);
         strcat(message," ");
         strcat(message,temp);
         strcat(message,"\n");
		 
   }

   if (((CReFoldDoc*) GetDocument())->ct.GetNumberofForbiddenPairs()>0) {
		strcat (message,"Prohibited Base Pairs=");
		for (it = 0; it<((CReFoldDoc*) GetDocument())->ct.GetNumberofForbiddenPairs();it++) {

 			itoa(((CReFoldDoc*) GetDocument())->ct.GetForbiddenPair5(it),temp,10);
			strcat(message,"  ");
			strcat(message,temp);
			strcat(message," - ");
			itoa(((CReFoldDoc*) GetDocument())->ct.GetForbiddenPair3(it),temp,10);
			strcat(message,temp);
			if (it%5==0) strcat(message,"\n");
			else strcat(message,", ");
		}
		strcat(message,"\n");
   }

   for (i=0;i<((CReFoldDoc*) GetDocument())->ct.nregion;i++) {

	    strcat(message,"Region from ");
		itoa(((CReFoldDoc*) GetDocument())->ct.start[i],temp,10);
		strcat(message,temp);
		strcat(message," to ");
		itoa(((CReFoldDoc*) GetDocument())->ct.stop[i],temp,10);
		strcat(message,temp);
		strcat(message,":\n");
		if (((CReFoldDoc*) GetDocument())->ct.rneighbors[i][0][0]>0) {
			strcat (message,"\tPaired neighbors:");
      
			for (it=0;it<100&&((CReFoldDoc*) GetDocument())->ct.rneighbors[i][it][0]>0;it++) {
				for (jt=0;jt<25&&((CReFoldDoc*) GetDocument())->ct.rneighbors[i][it][jt]>0;jt++) {
					strcat(message,tobase(((CReFoldDoc*) GetDocument())->ct.rneighbors[i][it][jt]));
				}
				strcat(message,"\n");
			 }
			 strcat(message,"\n");
		 
			}

		if (((CReFoldDoc*) GetDocument())->ct.rmin_gu[i]>0) {
			strcat (message,"\tMinimum # of G-U pairs=");
      
      		itoa(((CReFoldDoc*) GetDocument())->ct.rmin_gu[i],temp,10);
			strcat(message," ");
			strcat(message,temp);
			strcat(message,"\n");
		 
		}

		if (((CReFoldDoc*) GetDocument())->ct.rmin_g_or_u[i]>0) {
			strcat (message,"\tMinimum # of G's and U's in pairs=");
      
      		itoa(((CReFoldDoc*) GetDocument())->ct.rmin_g_or_u[i],temp,10);
			strcat(message," ");
			strcat(message,temp);
			strcat(message,"\n");
		 
		}

   }

   for (i=0;i<((CReFoldDoc*) GetDocument())->ct.nmicroarray;i++) {
		strcat(message,"Region from ");
		itoa(((CReFoldDoc*) GetDocument())->ct.microstart[i],temp,10);
		strcat(message,temp);
		strcat(message," to ");
		itoa(((CReFoldDoc*) GetDocument())->ct.microstop[i],temp,10);
		strcat(message,temp);
		strcat(message," has a minimum of ");
		itoa(((CReFoldDoc*) GetDocument())->ct.microunpair[i],temp,10);
		strcat(message,temp);
		strcat(message," unpaired nucs.\n");


   }

   

	if (((CReFoldDoc*) GetDocument())->ct.GetNumberofDoubles()==0&&((CReFoldDoc*) GetDocument())->ct.GetNumberofPairs()==0&&
		((CReFoldDoc*) GetDocument())->ct.GetNumberofSingles()==0&&((CReFoldDoc*) GetDocument())->ct.GetNumberofGU()==0&&
		((CReFoldDoc*) GetDocument())->ct.GetNumberofModified()==0&&((CReFoldDoc*) GetDocument())->ct.min_gu==0&&
		((CReFoldDoc*) GetDocument())->ct.min_g_or_u==0&&((CReFoldDoc*) GetDocument())->ct.nneighbors==0
		&&((CReFoldDoc*) GetDocument())->ct.GetNumberofForbiddenPairs()==0&&((CReFoldDoc*) GetDocument())->ct.nregion==0&&((CReFoldDoc*) GetDocument())->ct.nmicroarray==0) {

		strcpy(message,"There are no folding constraints.");
	}

   /*if (((CReFoldDoc*) GetDocument())->ct.ndbl>0) {
   	strcat (message,"Bases forced double:\n");
      for (it=1;it<=((CReFoldDoc*) GetDocument())->ct.ndbl;it++) {
      	itoa(((CReFoldDoc*) GetDocument())->ct.dbl[it],temp,10);
         strcat(message,"   ");
         strcat(message,temp);

         if (it%5==0) strcat(message,"\n");
		 else strcat(message,", ");
      }
	  strcat(message,"\n");
   }


   if (((CReFoldDoc*) GetDocument())->ct.npair>0) {
   	strcat(message,"Forced Base Pairs:\n");
   	for (it = 1; it<=((CReFoldDoc*) GetDocument())->ct.npair;it++) {

 			itoa(((CReFoldDoc*) GetDocument())->ct.pair[it][0],temp,10);
         strcat(message,"  ");
         strcat(message,temp);
         strcat(message," - ");
         itoa(((CReFoldDoc*) GetDocument())->ct.pair[it][1],temp,10);
         strcat(message,temp);
         if (it%5==0) strcat(message,"\n");
		 else strcat(message,", ");
      }
		strcat(message,"\n");
   }

   if (((CReFoldDoc*) GetDocument())->ct.nnopair>0) {
   	strcat (message,"Bases forced single:\n");
      for (it=1;it<=((CReFoldDoc*) GetDocument())->ct.nnopair;it++) {
      	itoa(((CReFoldDoc*) GetDocument())->ct.nopair[it],temp,10);
         strcat(message,"   ");
         strcat(message,temp);
         if (it%5==0) strcat(message,"\n");
		 else strcat(message,", ");
      }
	  strcat(message,"\n");
   }

   if (((CReFoldDoc*) GetDocument())->ct.ngu>0) {
   	strcat (message,"U's in GU pairs:\n");
      for (it=0;it<((CReFoldDoc*) GetDocument())->ct.ngu;it++) {
      	itoa(((CReFoldDoc*) GetDocument())->ct.gu[it],temp,10);
         strcat(message,"   ");
         strcat(message,temp);
         if (((it+1)%5==0)) strcat(message,"\n");
		 else strcat(message,", ");
      }
	  strcat(message,"\n");
   }

   if (((CReFoldDoc*) GetDocument())->ct.nmod>0) {
   	strcat (message,"Bases accessible to chemical modification:\n");
      for (it=1;it<=((CReFoldDoc*) GetDocument())->ct.nmod;it++) {
      	itoa(((CReFoldDoc*) GetDocument())->ct.mod[it],temp,10);
         strcat(message,"   ");
         strcat(message,temp);
         if (it%5==0) strcat(message,"\n");
		 else strcat(message,", ");
      }
	  strcat(message,"\n");
   }

   if (((CReFoldDoc*) GetDocument())->ct.min_gu>0) {
	  strcat (message,"Minimum # of G-U pairs=");
      
      	itoa(((CReFoldDoc*) GetDocument())->ct.min_gu,temp,10);
         strcat(message," ");
         strcat(message,temp);
         strcat(message,"\n");
		 
   }

   if (((CReFoldDoc*) GetDocument())->ct.neighbors[0][0]>0) {
	   strcat (message,"Paired neighbors:");
      
      for (it=0;it<100&&((CReFoldDoc*) GetDocument())->ct.neighbors[it][0]>0;it++) {
		  for (jt=0;jt<25&&((CReFoldDoc*) GetDocument())->ct.neighbors[it][jt]>0;jt++) {
			strcat(message,tobase(((CReFoldDoc*) GetDocument())->ct.neighbors[it][jt]));
		  }
		  strcat(message,"\n");
      }
	  strcat(message,"\n");
		 
   }

   if (((CReFoldDoc*) GetDocument())->ct.min_gu>0) {
	  strcat (message,"Minimum # of G-U pairs=");
      
      	itoa(((CReFoldDoc*) GetDocument())->ct.min_gu,temp,10);
         strcat(message," ");
         strcat(message,temp);
         strcat(message,"\n");
		 
   }

   if (((CReFoldDoc*) GetDocument())->ct.min_g_or_u>0) {
	  strcat (message,"Minimum # of G's and U's in pairs=");
      
      	itoa(((CReFoldDoc*) GetDocument())->ct.min_g_or_u,temp,10);
         strcat(message," ");
         strcat(message,temp);
         strcat(message,"\n");
		 
   }

   if (((CReFoldDoc*) GetDocument())->ct.nforbid>0) {
		strcat (message,"Prohibited Base Pairs=");
		for (it = 0; it<((CReFoldDoc*) GetDocument())->ct.nforbid;it++) {

 			itoa(((CReFoldDoc*) GetDocument())->ct.forbid[it][0],temp,10);
			strcat(message,"  ");
			strcat(message,temp);
			strcat(message," - ");
			itoa(((CReFoldDoc*) GetDocument())->ct.forbid[it][1],temp,10);
			strcat(message,temp);
			if (it%5==0) strcat(message,"\n");
			else strcat(message,", ");
		}
		strcat(message,"\n");
   }

   

	if (((CReFoldDoc*) GetDocument())->ct.ndbl==0&&((CReFoldDoc*) GetDocument())->ct.npair==0&&
		((CReFoldDoc*) GetDocument())->ct.nnopair==0&&((CReFoldDoc*) GetDocument())->ct.ngu==0&&
		((CReFoldDoc*) GetDocument())->ct.nmod==0&&((CReFoldDoc*) GetDocument())->ct.min_gu==0&&
		((CReFoldDoc*) GetDocument())->ct.min_g_or_u==0&&((CReFoldDoc*) GetDocument())->ct.neighbors[0][0]==0
		&&((CReFoldDoc*) GetDocument())->ct.nforbid==0) {

		strcpy(message,"There are no folding constraints.");
	}*/

   AfxMessageBox( message, 
			MB_OK|MB_ICONINFORMATION   );
	
}

void CRefoldView::OnForceSaveconstraints() 
{
		//The user is specifying the CT file name explicitly
	CFileDialog *filedialog;


	filedialog = new CFileDialog(FALSE,".con",NULL,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"Constraint Files (*.con)|*.con|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {


		outputconstraints(filedialog->GetPathName().GetBuffer(0),&(((CReFoldDoc*) GetDocument())->ct));

	

	}
	delete filedialog;
	
}
