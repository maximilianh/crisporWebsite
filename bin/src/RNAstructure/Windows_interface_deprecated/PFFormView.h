#pragma once

#include "folddoc.h"
#include "TProgressDialog.h"

UINT PFProc( LPVOID pParam );


//Foldobject:
struct CPFObject
{
	TProgressDialog *progress;
	structure *ct;
	datatable *data;

	double Temperature;
	
    CFormView *parent;
	
	char *savefile;
	bool subfold;
    

};

// CPFFormView form view

class CPFFormView : public CFormView
{
	DECLARE_DYNCREATE(CPFFormView)

protected:
	CPFFormView();           // protected constructor used by dynamic creation
	virtual ~CPFFormView();

public:
	enum { IDD = IDD_PFVIEW };
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedStart();
	afx_msg void OnBnClickedSequencebutton();
	afx_msg void OnBnClickedSavebutton();
	CFoldDoc *GetFoldDocument();
	CString sequencename;
	CString savename;
	void OnUpdate(CView*, LPARAM, CObject*);
	void OnInitialUpdate();
	CPFObject pfobject;
	TProgressDialog *progress;
	LRESULT DoneFolding(WPARAM wParam, LPARAM lParam);
	afx_msg void OnForceBasepair();
	afx_msg void OnForceChemicalmodification();
	afx_msg void OnForceDoublestranded();
	afx_msg void OnForceFmn();
	afx_msg void OnForceSingestranded();
	afx_msg void OnForceCurrent();
	afx_msg void OnForceReset();
	afx_msg void OnForceSaveconstraints();
	afx_msg void OnForceRestoreconstraints();
	afx_msg void OnForceForbidbasepairs();
	afx_msg void OnTemperature();
	afx_msg void OnForceMaximumpairingdistance();
	afx_msg void OnForceReadshapereactvity();
	afx_msg void OnReadShapeLinear();
	afx_msg void OnSubfold();

	bool subfold;
};


