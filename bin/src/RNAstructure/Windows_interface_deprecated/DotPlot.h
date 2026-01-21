#pragma once


// CDotPlot dialog

class CDotPlot : public CDialog
{
	DECLARE_DYNAMIC(CDotPlot)

public:
	CDotPlot(CWnd* pParent = NULL);   // standard constructor
	virtual ~CDotPlot();

// Dialog Data
	enum { IDD = IDD_DOTPLOT_DIALOG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedOk();
	afx_msg void OnBnClickedSavefile();
	afx_msg void OnBnClickedOutfile();
	CString savefilename;
	CString outfilename;
	float maxpercent;
	float maxenergy;
};
