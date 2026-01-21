#pragma once

#include "../src/structure.h"
#include "../src/algorithm.h"

// DynRefoldView form view

class DynRefoldView : public CFormView
{
	DECLARE_DYNCREATE(DynRefoldView)

protected:
	DynRefoldView();           // protected constructor used by dynamic creation
	virtual ~DynRefoldView();

public:
	enum { IDD = IDD_DYNREFOLDVIEW_FORM };
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedSavebutton();
	afx_msg void OnBnClickedCtbutton();
	afx_msg void OnBnClickedCtbutton2();
	afx_msg void OnBnClickedAlignment();
	afx_msg void OnBnClickedStart();
	int percent;
	int maxtracebacks;
	short window;
	short alignwindow;
	CString savefile;
	CString ct1file;
	CString ct2file;
	CString alignmentfile;
	bool alignallocated;
	short **align;
	int ct1size;
	virtual void OnInitialUpdate();
};


