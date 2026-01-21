// OligoObject.h: interface for the COligoObject class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_OLIGOOBJECT_H__251E7A92_BFC7_4D71_957D_1C452C4FF621__INCLUDED_)
#define AFX_OLIGOOBJECT_H__251E7A92_BFC7_4D71_957D_1C452C4FF621__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "../src/defines.h"
#include "../src/structure.h"


class COligoObject  
{
public:
	COligoObject(char *Datapath,char *Startpath,bool *Isdna, bool *Istargetdna, int *Option, 
		 int *Length,double *C,
		int *Usesub);
	bool AllocateTable();
	virtual ~COligoObject();

	char datapath[maxfil];
	char startpath[maxfil];
	bool isdna;
	bool istargetdna;
	int option;
	
	int length;
	double c;
	int usesub;
	int start;
	int stop;
	char outputfile[maxfil];
	int useprefilter;

	int **table,**numofsubstructures;
	bool allocated;

	//three filtering criteria for siRNA design:
	double asuf;
	double tofe;
	double fnnfe;
	float T;//The temperature of the calculation.

	int newtemp(bool isrna, datatable *freeenergy);



	structure ct;
	CFormView *parent;

	datatable data,*ddata,dhdata,*denthalpy;
	rddata *hybriddata;
    thermo *helixstack;

	TProgressDialog *update;


	//extra infrastructure for siRNA design:
	bool siRNA;//flag to indicate whether siRNA design is being performed
	bool *mask;



};

#endif // !defined(AFX_OLIGOOBJECT_H__251E7A92_BFC7_4D71_957D_1C452C4FF621__INCLUDED_)
