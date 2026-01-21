#pragma once


// CLinearSHAPEDialog dialog

class CLinearSHAPEDialog : public CDialog
{
	DECLARE_DYNAMIC(CLinearSHAPEDialog)

public:
	CLinearSHAPEDialog(structure *CT = NULL, CWnd* pParent = NULL);   // standard constructor
	virtual ~CLinearSHAPEDialog();
	structure *ct;

// Dialog Data
	enum { IDD = IDD_READSHAPEDIALOG_LINEAR };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	double Slope;
	double intercept;
	CString InputFilename;
	afx_msg void OnBnClickedOk();
	afx_msg void OnBnClickedCancel();
	afx_msg void OnBnClickedShapefilebutton();
};
