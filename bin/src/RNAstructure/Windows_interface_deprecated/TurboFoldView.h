#pragma once

#include "../TurboFold/TurboFold_object.h"
#include "TurboFoldDoc.h"
#include "afxwin.h"

//A function to kcik off a seperate thread for Multilign calculations
UINT TurboFoldProc( LPVOID pParam );

//MultilignObject, for passing information to a multilign thread:
struct CTurboFoldObject
{
	TProgressDialog *progress;//Need to add support for this.
    CFormView *parent;

	vector<string>  *savefilenames;
	vector<string>  *seqfilenames;
	vector<string>  *outputctfilenames;

	float Temperature;
	CString errormessage;
	char *datapath;

	int mode;//1 = MEA, 2= Threshold, 3= ProbKnot/TurboKnot

	double gamma;
	int iterations;

	float probability;//needs to be between 0 and 1

	double maxpercent, meagamma;//maxperecnt is between 0 and 100
	int maxstruct;
	int window;
	
	
	bool ISRNA;

	int knotiterations;
	int minlength;
	
	   

};

// CTurboFoldView form view

class CTurboFoldView : public CFormView
{
	DECLARE_DYNCREATE(CTurboFoldView)

protected:
	CTurboFoldView();           // protected constructor used by dynamic creation
	virtual ~CTurboFoldView();

public:
	enum { IDD = IDD_TURBOFOLD };
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
	afx_msg void OnSequencebutton();
	afx_msg void OnCtbutton();
	afx_msg void OnAdd();
	afx_msg void OnDelete();
	afx_msg void OnStart();
	CString m_sequence;
	CString m_ct;
	bool m_mea;
	bool m_pt;
	double m_gamma;
	int m_iterations;
	double m_percent;
	int m_maxstructures;
	int m_window;
	double m_meagamma;
	double m_threshold;
	int m_delete;
	bool started;

	float Temperature;
	CTurboFoldObject *object;


	CTurboFoldDoc *GetTurboFoldDocument(); 
	void OnInitialUpdate();
	void UpdateSequenceList();
	LRESULT Done(WPARAM wParam, LPARAM lParam);
	LRESULT DataError(WPARAM wParam, LPARAM lParam);
	LRESULT Error(WPARAM wParam, LPARAM lParam);

	


	TProgressDialog *progress;

	afx_msg void OnTemperature();
	CEdit C_SequenceList;
	afx_msg void OnEnChangeCtname();
	bool m_pk;
	// Number of iterations for a ProbKnot/TurboKnot assembly.
	int m_knotiterations;
	// Minimum helix length in a ProbKnot/TurboKnot prediction.
	int m_minlength;
};


