#pragma once


// CTemp_Dialog dialog

class CTemp_Dialog : public CDialog
{
	DECLARE_DYNCREATE(CTemp_Dialog)

public:
	CTemp_Dialog(float *temp = NULL, CWnd* pParent = NULL);   // standard constructor
	virtual ~CTemp_Dialog();
// Overrides
	void CTemp_Dialog::OnButtonOK();
	void CTemp_Dialog::OnButtonCancel();
	float *T,m_temp;



protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	virtual BOOL OnInitDialog();

	DECLARE_MESSAGE_MAP()
	
};
