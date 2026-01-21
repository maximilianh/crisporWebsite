#include "stdafx.h"
#include "boxplotview.h"
#include "plotdoc.h"
#include "../src/algorithm.h"
#include <math.h>
#include <iostream>
using namespace std;

IMPLEMENT_DYNCREATE(CBoxPlotView, CScrollView)

CBoxPlotView::CBoxPlotView(void)
{
}

CBoxPlotView::~CBoxPlotView(void)
{
}
BEGIN_MESSAGE_MAP(CBoxPlotView, CPlot)
	ON_COMMAND(ID_OUTPUT_OUTPUTPROBABLESTRUCTURE, OnOutputOutputprobablestructure)
	ON_COMMAND(ID_OUTPUT_OUTPUTTEXTFILE, OnOutputOutputtextfile)
	
END_MESSAGE_MAP()

void CBoxPlotView::OnOutputOutputprobablestructure()
{
	// output a single CT file with a structures with 99%, 97% , 95%, 90%, 80%, 70%, 60%, and 50% probable base pairs
	//The user is specifying the CT file name explicitly
	CFileDialog *filedialog;
	structure *ct;
	int count,i,j;
	double temp;
	
	
	


	filedialog = new CFileDialog(FALSE,".ct",0,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"CT Files (*.ct)|*.ct|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {
		//get the pointer to the document's ct file:
		ct = &(((PlotDoc *) GetDocument())->ct);

		//now construct the structures

		//zero the basepair arrays
		/*for (i=1;i<=ct->numofbases;i++) {
				ct->basepr[1][i] = 0;
				ct->basepr[2][i] = 0;
				ct->basepr[3][i] = 0;
				ct->basepr[4][i] = 0;
				ct->basepr[5][i] = 0;
				ct->basepr[6][i] = 0;
				ct->basepr[7][i] = 0;
				ct->basepr[8][i] = 0;
				
			
		}*/

		//Allocate the needed structures:
		for (count=1;count<=8;count++) {
			ct->AddStructure();
		}

		for (i=1;i<ct->GetSequenceLength();i++) {
			for (j=i+1;j<=ct->GetSequenceLength();j++) {

				//double temp2 = ((PlotDoc *)GetDocument())->array[j][i];
				
				if (((PlotDoc *)GetDocument())->arrayvalues[j][i]>0) {  //a value > 0 is stored for pairs that cannot occur
				for (count=1;count<=8;count++) {

					//temp = 	pow(10,-((PlotDoc *)GetDocument())->array[j][i]);

					if (count==1) {
						if (pow(10,-((PlotDoc *)GetDocument())->arrayvalues[j][i])>=.99) {
						//ct->basepr[count][i]=j;
						//ct->basepr[count][j]=i;
						ct->SetPair(i,j,count);
						}
					}
					else if (count==2) {
							if (pow(10,-((PlotDoc *)GetDocument())->arrayvalues[j][i])>=.97) {
						//ct->basepr[count][i]=j;
						//ct->basepr[count][j]=i;
								ct->SetPair(i,j,count);
							}
					}
					else if (count==3) {
						if (pow(10,-((PlotDoc *)GetDocument())->arrayvalues[j][i])>=.95) {
						//ct->basepr[count][i]=j;
						//ct->basepr[count][j]=i;
							ct->SetPair(i,j,count);
						}
					}
					else if (count==4) {
						if (pow(10,-((PlotDoc *)GetDocument())->arrayvalues[j][i])>=.90) {
						//ct->basepr[count][i]=j;
						//ct->basepr[count][j]=i;
							ct->SetPair(i,j,count);
						}
					}
					else if (count==5) {
						if (pow(10,-((PlotDoc *)GetDocument())->arrayvalues[j][i])>=.80) {
						//ct->basepr[count][i]=j;
						//ct->basepr[count][j]=i;
							ct->SetPair(i,j,count);
						}
					}
					else if (count==6) {
						if (pow(10,-((PlotDoc *)GetDocument())->arrayvalues[j][i])>=.70) {
						//ct->basepr[count][i]=j;
						//ct->basepr[count][j]=i;
							ct->SetPair(i,j,count);
						}
					}
					else if (count==7) {
						if (pow(10,-((PlotDoc *)GetDocument())->arrayvalues[j][i])>=.60) {
						//ct->basepr[count][i]=j;
						//ct->basepr[count][j]=i;
							ct->SetPair(i,j,count);
						}
					}
					else if (count==8) {
						if (pow(10,-((PlotDoc *)GetDocument())->arrayvalues[j][i])>.50) {
						//ct->basepr[count][i]=j;
						//ct->basepr[count][j]=i;
							ct->SetPair(i,j,count);
						}
					}
					
				}
				}
			}
			
		}

		
		string ctlabel;
		ctlabel = " >=99% probable pairs " + ct->GetSequenceLabel() ;
		ct->SetCtLabel(ctlabel,1);

		ctlabel = " >=97% probable pairs " + ct->GetSequenceLabel() ;
		ct->SetCtLabel(ctlabel,2);

		ctlabel = " >=95% probable pairs " + ct->GetSequenceLabel()  ;
		ct->SetCtLabel(ctlabel,3);

		ctlabel = " >=90% probable pairs " + ct->GetSequenceLabel()  ;
		ct->SetCtLabel(ctlabel,4);

		ctlabel = " >=80% probable pairs " + ct->GetSequenceLabel()  ;
		ct->SetCtLabel(ctlabel,5);

		ctlabel = " >=70% probable pairs " + ct->GetSequenceLabel()  ;
		ct->SetCtLabel(ctlabel,6);

		ctlabel = " >=60% probable pairs " + ct->GetSequenceLabel()  ;
		ct->SetCtLabel(ctlabel,7);

		ctlabel = " >50% probable pairs " + ct->GetSequenceLabel()  ;
		ct->SetCtLabel(ctlabel,8);



		/*strcpy(ct->ctlabel[9],ct->ctlabel[1]);
		strcpy(ct->ctlabel[2]," >=97% probable pairs ");
		strcpy(ct->ctlabel[3]," >=95% probable pairs ");
		strcpy(ct->ctlabel[4]," >=90% probable pairs ");
		strcpy(ct->ctlabel[5]," >=80% probable pairs ");
		strcpy(ct->ctlabel[6]," >=70% probable pairs ");
		strcpy(ct->ctlabel[7]," >=60% probable pairs ");
		strcpy(ct->ctlabel[8]," >50% probable pairs ");
		strcat(ct->ctlabel[2],ct->ctlabel[1]);
		strcat(ct->ctlabel[3],ct->ctlabel[1]);
		strcat(ct->ctlabel[4],ct->ctlabel[1]);
		strcat(ct->ctlabel[5],ct->ctlabel[1]);
		strcat(ct->ctlabel[6],ct->ctlabel[1]);
		strcat(ct->ctlabel[7],ct->ctlabel[1]);
		strcat(ct->ctlabel[8],ct->ctlabel[1]);
		strcpy(ct->ctlabel[1]," >=99% probable pairs ");
		strcat(ct->ctlabel[1],ct->ctlabel[9]);*/

		//ct->numofstructures = 8;
		ct->ctout((filedialog->GetPathName()).GetBuffer(0));
				

	}

	//offer to display the predicted structures if this is not subfolding
	CDialog *finished = new CDialog(IDD_FINISHED,this);

	if(finished->DoModal()==IDOK) ( ((PlotDoc *)GetDocument())->app)->Draw((filedialog->GetPathName()).GetBuffer(0));
	

	delete finished;

	delete filedialog;



}

void CBoxPlotView::OnOutputOutputtextfile()
{
	//Output a Tab-delimitted text file with all base pair probabilities

	ofstream out;
	CFileDialog *filedialog;
	structure *ct;
	int i,j;

	//Get a file name first

	filedialog = new CFileDialog(FALSE,".ct",0,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"Dot Plot Files (*.dp)|*.dp|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {
		ct = &(((PlotDoc *) GetDocument())->ct);
		out.open((filedialog->GetPathName()).GetBuffer(0));
		out <<ct->GetSequenceLength()<<"\n";
		out << "i" << "\t"<<"j"<<"\t"<<"-log10(Probability)"<<"\n";
		for (i=1;i<ct->GetSequenceLength();i++) {
			for (j=i+minloop+1;j<=ct->GetSequenceLength();j++) {
				if (plotdoc->arrayvalues[j][i]<=plotdoc->colorrange[plotdoc->colors]&&plotdoc->arrayvalues[j][i]>=plotdoc->colorrange[0])
					out << i << "\t"<<j<<"\t"<<(((PlotDoc *)GetDocument())->arrayvalues[j][i])<<"\n"; 

			}
		}
		out.close();
	}

}

