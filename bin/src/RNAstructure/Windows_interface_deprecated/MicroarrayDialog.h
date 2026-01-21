#pragma once
#include "../src/structure.h"

// CMicroarrayDialog dialog

class CMicroarrayDialog : public CDialog
{
	DECLARE_DYNAMIC(CMicroarrayDialog)

public:
	CMicroarrayDialog(CWnd* pParent = NULL);   // standard constructor
	CMicroarrayDialog(structure *CT, CWnd* pParent = NULL);
	virtual ~CMicroarrayDialog();
	structure *ct;

// Dialog Data
	enum { IDD = IDD_MICROARRAYDIALOG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedOk();
	int m_start;
	int m_stop;
	int m_min;
};
