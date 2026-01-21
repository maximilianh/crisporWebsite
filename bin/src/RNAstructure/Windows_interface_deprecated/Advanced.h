#if !defined(AFX_ADVANCED_H__8F296020_965F_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_ADVANCED_H__8F296020_965F_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// Advanced.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CAdvanced dialog

class CAdvanced : public CDialog
{
// Construction
public:
	CAdvanced(short *dynnumber, short *dynpercent, short *dynwindow,CWnd* pParent = NULL);   // standard constructor
	short *dynnumber,*dynpercent,*dynwindow;

// Dialog Data
	//{{AFX_DATA(CAdvanced)
	enum { IDD = IDD_ADVANCED };
	short	*m_dnumber;
	short	*m_dwindow;
	short	*m_dpercent;
	//}}AFX_DATA
	short savenumber,savewindow,savepercent;

	void OnChangeParameter();
	void OnCancel();


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CAdvanced)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CAdvanced)
	virtual void OnOK();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_ADVANCED_H__8F296020_965F_11D4_9F32_00C0F02A5F5D__INCLUDED_)
