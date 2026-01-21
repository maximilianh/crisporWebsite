#pragma once

#include "TProgressDialog.h"
#include "../src/phmm/utils/xmath/log/xlog_math.h"

//PKobject:
struct CPKObject
{
	TProgressDialog *progress;
	
	
	
	int minlength;
	int iterations;
	
	
    CFormView *parent;
	char *ctoutfile;
	char *savefile;
	
   

};



// CProbKnotView form view

class CProbKnotView : public CFormView
{
	DECLARE_DYNCREATE(CProbKnotView)

protected:
	CProbKnotView();           // protected constructor used by dynamic creation
	virtual ~CProbKnotView();

public:
	enum { IDD = IDD_PROBKNOTVIEW };
	void OnUpdate(CView*, LPARAM, CObject*);
#ifdef _DEBUG
	virtual void AssertValid() const;
#ifndef _WIN32_WCE
	virtual void Dump(CDumpContext& dc) const;
#endif
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	CString m_savefile;
public:
	CString m_ctfilename;
public:
	afx_msg void OnBnClickedOpensavebutton();
public:
	afx_msg void OnBnClickedCtbutton();
public:
	afx_msg void OnBnClickedStart();


private:
	CFoldDoc *GetFoldDocument();
	bool started;
public:
	int m_iterations;
	int m_ml;
	LRESULT DoneFolding(WPARAM wParam, LPARAM lParam);
	TProgressDialog *progress;
	CPKObject pkobject;
};


