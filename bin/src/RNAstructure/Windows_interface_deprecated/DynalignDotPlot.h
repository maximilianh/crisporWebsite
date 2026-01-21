#pragma once


// CDynalignDotPlot dialog

class CDynalignDotPlot : public CDialog
{
	DECLARE_DYNAMIC(CDynalignDotPlot)

public:
	CDynalignDotPlot(CWnd* pParent = NULL);   // standard constructor
	virtual ~CDynalignDotPlot();

// Dialog Data
	enum { IDD = IDD_DYNALIGNDOTPLOT_DIALOG };

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
