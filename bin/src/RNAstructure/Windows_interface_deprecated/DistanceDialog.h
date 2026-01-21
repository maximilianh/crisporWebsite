#pragma once

#include "../src/structure.h"


// CDistanceDialog dialog

class CDistanceDialog : public CDialog
{
	DECLARE_DYNAMIC(CDistanceDialog)

public:
	CDistanceDialog(structure *CT = NULL, CWnd* pParent = NULL);   // standard constructor
	virtual ~CDistanceDialog();
	structure *ct;
	

// Dialog Data
	enum { IDD = IDD_DISTANCEDIALOG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
	BOOL OnInitDialog();
public:
	int m_distance;
	
	afx_msg void OnBnClickedOk();
};
