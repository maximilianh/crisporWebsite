#pragma once


// CAlignForceDialog dialog

class CAlignForceDialog : public CDialog
{
	DECLARE_DYNAMIC(CAlignForceDialog)

public:
	CAlignForceDialog(short **Alignforce= NULL, int Length1=0, int Length2=0, CWnd* pParent = NULL);   // standard constructor
	virtual ~CAlignForceDialog();
	short **alignforce;
	int length1,length2; //The lengths of the sequence1 and sequence2, respectively, passed to constructor

// Dialog Data
	enum { IDD = IDD_ALIGNDIALOG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	int m_ns1;
	int m_ns2;
	afx_msg void OnBnClickedOkandopen();
	afx_msg void OnBnClickedOk();
};
