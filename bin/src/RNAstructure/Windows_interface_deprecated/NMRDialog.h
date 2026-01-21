#pragma once


// NMRDialog dialog

class NMRDialog : public CDialog
{
	DECLARE_DYNAMIC(NMRDialog)

public:
	NMRDialog(CWnd* pParent = NULL);   // standard constructor
	NMRDialog(structure *CT, CWnd* pParent = NULL); //constructor to be used
	virtual ~NMRDialog();
	structure *ct;

// Dialog Data
	enum { IDD = IDD_NMR };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedOk();
	afx_msg void OnBnClickedCancel();
	int min_g_or_u;
	int min_gu_pair;
	CString neighbors;
};
