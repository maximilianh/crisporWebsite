#pragma once


// CInternalDialog dialog

class CInternalDialog : public CDialog
{
	DECLARE_DYNAMIC(CInternalDialog)

public:
	CInternalDialog(int *Maxinternal=NULL, CWnd* pParent = NULL);   // standard constructor
	virtual ~CInternalDialog();
	int *maxinternal,m_max;
	void OnButtonOK();
	void OnButtonCancel();

// Dialog Data
	enum { IDD = IDD_INTERNALDIALOG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
};
