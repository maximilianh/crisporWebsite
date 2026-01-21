#pragma once


// CRegionalNMR dialog

class CRegionalNMR : public CDialog
{
	DECLARE_DYNAMIC(CRegionalNMR)

public:
	CRegionalNMR(CWnd* pParent = NULL);   // standard constructor
	CRegionalNMR(structure *CT, CWnd* pParent = NULL); //constructor to be used
	virtual ~CRegionalNMR();
	structure *ct;

// Dialog Data
	enum { IDD = IDD_REGIONAL_NMR };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	int min_g_or_u;
	int min_gu_pair;
	CString neighbors;
	int start,stop;
	afx_msg void OnBnClickedOk();
	afx_msg void OnBnClickedCancel();
};
