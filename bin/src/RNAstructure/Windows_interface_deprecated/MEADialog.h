#pragma once

#include "TProgressDialog.h"

//MEAobject:
struct CMEAObject
{
	TProgressDialog *progress;
	structure *ct;
	
	double percent;
	int structures;
	double gamma;
	int window;
	
	
    CFormView *parent;
	char *ctoutfile;
	char *savefile;
	
   

};


// CMEADialog dialog

class CMEADialog : public CFormView
{
	DECLARE_DYNCREATE(CMEADialog)

public:
	CMEADialog(CWnd* pParent = NULL);   // standard constructor
	LRESULT DoneFolding(WPARAM wParam, LPARAM lParam);
	CFoldDoc *GetFoldDocument();
	void OnUpdate(CView*, LPARAM, CObject*);
	virtual ~CMEADialog();
	#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

// Dialog Data
	enum { IDD = IDD_MEA };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedPartition();

	afx_msg void OnBnClickedCtbutton();

	afx_msg void OnBnClickedStart();

	CString m_SEQUENCE;

	CString m_CT;

private:
	bool started;
	CMEAObject meaobject;//Use a stochastic object to pass information to a seperate thread
	double m_percent;
	int m_structures;
	int m_window;
	double m_gamma;
	TProgressDialog *progress;
};
