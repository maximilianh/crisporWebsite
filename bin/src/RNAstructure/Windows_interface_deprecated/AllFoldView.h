#pragma once

#include "../src/algorithm.h"
#include "FoldDoc.h"
#include "TProgressDialog.h"

UINT AllFoldProc( LPVOID pParam );


//Foldobject:
struct CAllFoldObject
{
	
	TProgressDialog *progress;
	structure *ct;
	datatable *data;
	int m_percent;
	int m_abs;
    CFormView *parent;
	char *ctoutfile;
	char *savefile;
	
    

};


// AllFoldView dialog

class AllFoldView : public CFormView
{
	DECLARE_DYNCREATE(AllFoldView)

public:
	AllFoldView(CWnd* pParent = NULL);   // standard constructor
	virtual ~AllFoldView();

// Dialog Data
	enum { IDD = IDD_ALLFOLDS };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedSequencebutton();
	afx_msg void OnBnClickedCtbutton();
	void NewSequence();
	CFoldDoc *GetFoldDocument();
	afx_msg void OnBnClickedStart();
	BOOL m_save;
	short m_percent;
	afx_msg void OnForceBasepair();
	afx_msg void OnForceChemicalmodification();
	afx_msg void OnForceDoublestranded();
	afx_msg void OnForceFmn();
	afx_msg void OnForceNmrconstraints();
	afx_msg void OnForceSingestranded();
	afx_msg void OnForceForbidbasepairs169();
	afx_msg void OnForceCurrent();
	afx_msg void OnForceReset();
	afx_msg void OnForceSaveconstraints();
	afx_msg void OnForceRestoreconstraints();
	CString m_sequencename;
	CString m_ctname;
	float m_abs;
	void OnInitialUpdate();
	TProgressDialog *progress;
	CAllFoldObject foldobject;
	LRESULT DoneFolding(WPARAM wParam, LPARAM lParam);
	afx_msg void OnForceRegionalnmrconstraints();
	afx_msg void OnForceMicroarrayconstraints();
	void OnTemperature();
};
