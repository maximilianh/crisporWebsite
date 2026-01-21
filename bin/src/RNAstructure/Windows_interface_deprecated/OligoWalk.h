#if !defined(AFX_OLIGOWALK_H__782B7593_077C_4035_825C_A37460F2A543__INCLUDED_)
#define AFX_OLIGOWALK_H__782B7593_077C_4035_825C_A37460F2A543__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// OligoWalk.h : header file
//

#include "OligoDoc.h"

/////////////////////////////////////////////////////////////////////////////
// COligoWalk dialog

class COligoWalk : public CFormView
{
// Construction
protected:
	COligoWalk(UINT i=COligoWalk::IDD);   // standard constructor
	DECLARE_DYNCREATE(COligoWalk)
public:
	
	bool started;
	void OnInitialUpdate();
	
	void OnLengthSpin(NMHDR* pNMHDR, LRESULT* pResult);
	void OnStartSpin(NMHDR* pNMHDR, LRESULT* pResult);
	void OnStopSpin(NMHDR* pNMHDR, LRESULT* pResult);
	bool ctloaded;
	void OnTemperature();//User is changing the temperature of calculation.
	
	//These three parameters are used for siRNA design
	double m_asuf,m_tofe,m_fnnfe;

	COligoDoc* GetOligoDocument();

	LRESULT DisplayOligoWalk(WPARAM wParam, LPARAM lParam);

	CSpinButtonCtrl m_LengthSpin,m_StartSpin,m_StopSpin;

// Dialog Data
	//{{AFX_DATA(COligoWalk)
	enum { IDD = IDD_OLIGOWALK };
	int		m_length;
	int		m_start;
	int		m_stop;
	CString	m_ctname;
	CString	m_reportname;
	double	m_concentration;
	//}}AFX_DATA

	// Generated message map functions
	//{{AFX_MSG(COligoWalk)
	afx_msg void OnCtbutton();
	afx_msg void OnReportbutton();
	virtual void OnStart();
	afx_msg void OnChangeLength();
	afx_msg void OnChangeStart();
	afx_msg void OnChangeStop();
	//}}AFX_MSG

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(COligoWalk)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_OLIGOWALK_H__782B7593_077C_4035_825C_A37460F2A543__INCLUDED_)
