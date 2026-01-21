#if !defined(AFX_DYNALIGN_H__44DAAC49_F63B_4B50_AB83_8415247E630E__INCLUDED_)
#define AFX_DYNALIGN_H__44DAAC49_F63B_4B50_AB83_8415247E630E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// Dynalign.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CDynalign form view

#ifndef __AFXEXT_H__
#include <afxext.h>
#endif


#include "dyndoc.h"




//Dynobject:
struct CDynObject
{
	TProgressDialog *progress;
	structure *ct1;
	structure *ct2;
	datatable *data;
    CFormView *parent;
	char *ctoutfile;
	char *ctoutfile2;
	char *aoutfile;
	bool singleinsert;
	short m,dggap;
	short structwindow, alignwindow, maxtracebacks,percent;
	char *savefile;
	short **forcealign;
	char *datapath;
    

};

class CDynalign : public CFormView
{
protected:
	CDynalign();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CDynalign)

// Form Data
public:
	//{{AFX_DATA(CDynalign)
	enum { IDD = IDD_DYNALIGN_FORM };
	CString	m_ctname;
	CString	m_ctname2;
	CString	m_sequence2name;
	CString	m_sequencename;
	short	m_M;
	CString	m_alignment;
	BOOL	m_singleinsert;
	double	m_dggap;
	//}}AFX_DATA

	void OnInitialUpdate();
	CDynDoc *GetDynDocument();
	void NewSequences();
	TProgressDialog *progress;
	CDynObject dynobject;
	LRESULT DynDone(WPARAM wParam, LPARAM lParam);
	short **forcealign;
	void OnTemperature();

// Attributes
public:

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CDynalign)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	virtual ~CDynalign();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
	//{{AFX_MSG(CDynalign)
	afx_msg void OnSequencebutton();
	afx_msg void OnSequence2button();
	afx_msg void OnCtbutton();
	afx_msg void OnCtbutton2();
	afx_msg void OnAlignmentbutton();
	afx_msg void OnStart();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:
	// Maximum percent change in free energy when generating suboptimal structures.
	short percent;
	// Maximum number of suboptimal structures.
	short maxtracebacks;
	short structwindow;
	short alignwindow;
	BOOL savefilecheck;
	afx_msg void OnS1Bp();
	afx_msg void OnS1Cm();
	afx_msg void OnS1Ds();
	afx_msg void OnS1Fmn();
	afx_msg void OnS1Ss();
	afx_msg void OnS1Pb();
	afx_msg void OnS1Current();
	afx_msg void OnS1Reset();
	afx_msg void OnS1Save();
	afx_msg void OnS1Restore();
	afx_msg void OnS2Bp();
	afx_msg void OnS2Cm();
	afx_msg void OnS2Ds();
	afx_msg void OnS2Fmn();
	afx_msg void OnS2Ss();
	afx_msg void OnS2Pb();
	afx_msg void OnS2Current();
	afx_msg void OnS2Reset();
	afx_msg void OnS2Save();
	afx_msg void OnS2Restore();
	afx_msg void OnConstraintsforalignmentForcealignment();
	afx_msg void OnConstraintsforalignmentShowcurrentalignmentconstraints();
	afx_msg void OnConstraintsforalignmentResetalignmentcontraints();
	afx_msg void OnConstraintsforalignmentSavealignment();
	afx_msg void OnConstraintsforalignmentRestorealignment();
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_DYNALIGN_H__44DAAC49_F63B_4B50_AB83_8415247E630E__INCLUDED_)
