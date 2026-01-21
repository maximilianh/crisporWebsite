#pragma once


// CBoxPlotOpen dialog

class CBoxPlotOpen : public CDialog
{
	DECLARE_DYNAMIC(CBoxPlotOpen)

public:
	CBoxPlotOpen(CWnd* pParent = NULL);   // standard constructor
	CBoxPlotOpen(CString *Filename, CWnd* pParent = NULL);
	virtual ~CBoxPlotOpen();
	CString *filename;

// Dialog Data
	enum { IDD = IDD_BOXPLOTOPEN };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedOpensave();
	afx_msg void OnBnClickedOk();
	CString pfsave;
};
