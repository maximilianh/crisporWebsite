#pragma once
#include "../src/structure.h"
UINT StochasticProc( LPVOID pParam );
#include "TProgressDialog.h"

//Stochasticobject:
struct CStochasticObject
{
	TProgressDialog *progress;
	structure *ct;
	datatable *data;
	int m_seed;
	int m_structures;
	
    CFormView *parent;
	char *ctoutfile;
	char *savefile;
	
   

};

// CStochasticView dialog

class CStochasticView : public CFormView
{
	DECLARE_DYNCREATE(CStochasticView)

public:
	CStochasticView(CWnd* pParent = NULL);   // standard constructor
	virtual ~CStochasticView();
	CFoldDoc *GetFoldDocument();
	LRESULT DoneFolding(WPARAM wParam, LPARAM lParam);
	bool started;
	void OnUpdate(CView*, LPARAM, CObject*);


// Dialog Data
	enum { IDD = IDD_STOCHASTIC };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedPartition();
	afx_msg void OnBnClickedCtbutton();
	int m_structures;
	int m_seed;
	afx_msg void OnBnClickedStart();
	CString m_pfs;
	CString m_ct;
	CStochasticObject stochasticobject;
	TProgressDialog *progress;
};
