#pragma once

#include "../src/structure.h"


// CReadSHAPE dialog

class CReadSHAPE : public CDialog
{
	DECLARE_DYNAMIC(CReadSHAPE)

public:
	CReadSHAPE(structure *CT = NULL, CWnd* pParent = NULL);   // standard constructor
	virtual ~CReadSHAPE();
	structure *ct;

// Dialog Data
	enum { IDD = IDD_READSHAPEDIALOG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	CString InputFilename;
	float ForceSingleStranded;
	float ChemicalModification;
	afx_msg void OnBnClickedOk();
	afx_msg void OnBnClickedShapefilebutton();
};
