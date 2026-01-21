#pragma once
#include "afxwin.h"
#include "MultilignDoc.h"
#include "../RNA_class/Multilign_object.h"


//A function to kcik off a seperate thread for Multilign calculations
UINT MultilignProc( LPVOID pParam );

//MultilignObject, for passing information to a multilign thread:
struct CMultilignObject
{
	TProgressDialog *progress;
    CFormView *parent;
	vector<vector<string> > *filenames;
	char *alignment;
	double maxpercent;
	int maxstruct;
	int structurewindow;
	int alignmentwindow;
	double gap;
	bool singlebp;
	bool generatesave;
	int iterations;
	int maxpairs;
	float maxdsvchange;
	float Temperature;
	char *datapath;
	bool ISRNA;
	CString errormessage;
	   

};



// MultilignView form view

class MultilignView : public CFormView
{
	DECLARE_DYNCREATE(MultilignView)

protected:
	MultilignView();           // protected constructor used by dynamic creation
	virtual ~MultilignView();

public:
	enum { IDD = IDD_MULTILIGN };
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
	afx_msg void OnBnClickedSequencebutton();
	afx_msg void OnBnClickedCtbutton();
	CString m_sequence;
	CString m_ct;
	afx_msg void OnBnClickedAdd2();
	afx_msg void OnBnClickedAlignmentbutton();
	CString m_alignment;
	double m_maxpercent;
	int m_maxstructures;
	int m_structurewindow;
	int m_alignmentwindow;
	double m_gap;
	BOOL m_singlebp;
	BOOL m_generatesave;
	int m_iterations;
	int m_maxpairs;
	float m_maxdsvchange;
	afx_msg void OnBnClickedStart();
	int m_deletenumber;
	afx_msg void OnBnClickedDelete();
	CEdit C_SequenceList;
	float Temperature;
	CMultilignObject *object;
	bool started;


	CMultilignDoc *GetMultilignDocument(); 
	void OnInitialUpdate();
	void UpdateSequenceList();
	LRESULT Done(WPARAM wParam, LPARAM lParam);
	LRESULT DataError(WPARAM wParam, LPARAM lParam);
	LRESULT Error(WPARAM wParam, LPARAM lParam);


	TProgressDialog *progress;
public:
	afx_msg void OnTemperature();
};


