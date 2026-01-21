#if !defined(AFX_FOLDVIEW_H__941BF522_9231_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_FOLDVIEW_H__941BF522_9231_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// FoldView.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CFoldView form view

#ifndef __AFXEXT_H__
#include <afxext.h>
#endif


#include "../src/algorithm.h"
#include "FoldDoc.h"


UINT FoldProc( LPVOID pParam );
void currentconstraints(structure *ct);
void resetconstraints(structure *ct);


//Foldobject:
struct CFoldObject
{
	TProgressDialog *progress;
	structure *ct;
	datatable *data;
	int m_number;
	int dynpercent;
	int m_percent;
	int dynwindow;
	int m_window;
	int dynnumber;
    CFormView *parent;
	char *ctoutfile;
	char *savefile;
	int maximuminternalloopsize;
	bool subfold;
	
    

};


class CFoldView : public CFormView
{
protected:
	CFoldView();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CFoldView)

// Form Data
public:
	//{{AFX_DATA(CFoldView)
	enum { IDD = IDD_FOLDVIEW_FORM };
	CString	m_ctname;
	short	m_number;
	short	m_percent;
	BOOL	m_save;
	short	m_window;
	CString	m_sequencename;
	//}}AFX_DATA

	short dynwindow;
	short dynnumber;
	short dynpercent;
	TProgressDialog *progress;
	CFoldObject foldobject;
	void OnInitialUpdate();
	bool started;
	void OnTemperature();
	bool subfold;  //if true, the fold all subfragments from the 5' end.


// Attributes
public:
	void OnUpdate(CView*, LPARAM, CObject*);
	LRESULT DoneFolding(WPARAM wParam, LPARAM lParam); 
	CFoldDoc *GetFoldDocument();
// Operations
public:

private:
	void NewSequence();
	
	

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CFoldView)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	virtual ~CFoldView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
	//{{AFX_MSG(CFoldView)
	afx_msg void OnSequencebutton();
	afx_msg void OnCtbutton();
	afx_msg void OnStart();
	afx_msg void OnForceBasepair();
	afx_msg void OnForceCurrent();
	afx_msg void OnForceDoublestranded();
	afx_msg void OnForceReset();
	afx_msg void OnForceRestoreconstraints();
	afx_msg void OnForceSaveconstraints();
	afx_msg void OnForceSingestranded();
	//afx_msg void OnAdvancedSuboptimalstructures();
	afx_msg void OnForceFmn();
	afx_msg void OnChangeParameter();
	afx_msg void OnForceChemical();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnForceNmrconstraints();
	afx_msg void OnForceForbidbasepairs169();
	afx_msg void OnMaximumloop();
	afx_msg void OnForceReadshapereactvity();
	afx_msg void OnForceMaximumpairingdistance();
	afx_msg void OnReadShapeLinear();
	afx_msg void OnSubfold();
	
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_FOLDVIEW_H__941BF522_9231_11D4_9F32_00C0F02A5F5D__INCLUDED_)
