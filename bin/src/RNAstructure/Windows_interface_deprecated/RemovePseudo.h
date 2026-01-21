#pragma once

#include "FoldDoc.h"
#include "TProgressDialog.h"
#include "AllFoldView.h" 

UINT BPProc( LPVOID pParam );


//Foldobject:
struct CBPObject
{
	
	TProgressDialog *progress;
	structure *ct;
	datatable *data;
    CFormView *parent;
	char *ctoutfile;
	
	
    

};

// CRemovePseudo form view

class CRemovePseudo : public CFormView
{
	DECLARE_DYNCREATE(CRemovePseudo)

protected:
	CRemovePseudo();           // protected constructor used by dynamic creation
	virtual ~CRemovePseudo();

public:
	enum { IDD = IDD_REMOVEPSEUDO };
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedCtin();
	afx_msg void OnBnClickedCtbutton();
	afx_msg void OnBnClickedStart();

	CString m_ctinname;
	CString m_ctoutname;
	CFoldDoc *GetFoldDocument();
	TProgressDialog *progress;
	CAllFoldObject foldobject;
	LRESULT DoneFolding(WPARAM wParam, LPARAM lParam);
	void OnInitialUpdate();
	void OnTemperature();
};


