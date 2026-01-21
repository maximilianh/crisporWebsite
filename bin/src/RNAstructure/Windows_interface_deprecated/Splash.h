#if !defined(AFX_SPLASH_H__58413FA6_E53A_4FB8_8BA6_EA57A77F1BD8__INCLUDED_)
#define AFX_SPLASH_H__58413FA6_E53A_4FB8_8BA6_EA57A77F1BD8__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// Splash.h : header file
//
struct CTimeObject {

	CDialog *splash;

};
/////////////////////////////////////////////////////////////////////////////
// CSplash dialog

class CSplash : public CDialog
{
// Construction
public:
	CSplash(CWnd* pParent = NULL);   // standard constructor
	~CSplash();
	BOOL OnInitDialog( );
	LRESULT Done(WPARAM wParam, LPARAM lParam);
	CTimeObject time;
	CWnd* parent;
	DECLARE_DYNAMIC(CSplash)
	void OnOK();

// Dialog Data
	//{{AFX_DATA(CSplash)
	enum { IDD = IDD_SPLASH };
		// NOTE: the ClassWizard will add data members here
	//}}AFX_DATA


// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CSplash)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	virtual void PostNcDestroy();
	//}}AFX_VIRTUAL

// Implementation
protected:

	// Generated message map functions
	//{{AFX_MSG(CSplash)
	afx_msg void OnLButtonDblClk(UINT nFlags, CPoint point);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);

	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_SPLASH_H__58413FA6_E53A_4FB8_8BA6_EA57A77F1BD8__INCLUDED_)
