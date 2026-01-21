#pragma once

#include "oligoscreendoc.h"

//FOligoScreenObject to pass all the data needed in a seperate thread:
struct COligoScreenObject
{
	TProgressDialog *progress;
	char *infilename;
	char *outfilename;
	bool isRNA;
	char *path;
	CFormView *parent;
	float T;
	

};

// COligoScreen form view

class COligoScreen : public CFormView
{
	DECLARE_DYNCREATE(COligoScreen)

protected:
	COligoScreen();           // protected constructor used by dynamic creation
	virtual ~COligoScreen();

public:
	enum { IDD = IDD_OligoScreen };
	COligoScreenDoc* GetOligoScreenDocument();
	COligoScreenObject oligoscreenobject;
	TProgressDialog *progress;
	LRESULT DoneFolding(WPARAM wParam, LPARAM lParam); 
	void OnTemperature();
	float T;//Tempaerature of calculation.
	//char path[_MAX_PATH];
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	CString m_listname;
	CString m_outname;
	void OnInitialUpdate(void);
	afx_msg void OnBnClickedOk();
	afx_msg void OnBnClickedList();
	afx_msg void OnBnClickedOut();
};


