#pragma once
#include "PFFormView.h"


// CBimolPFView form view

class CBimolPFView : public CFormView
{
	DECLARE_DYNCREATE(CBimolPFView)

protected:
	CBimolPFView();           // protected constructor used by dynamic creation
	virtual ~CBimolPFView();

public:
	enum { IDD = IDD_BIMOL_PF };
	CFoldDoc *GetFoldDocument();
	CPFObject pfobject;
	TProgressDialog *progress;
	LRESULT DoneFolding(WPARAM wParam, LPARAM lParam);
	void OnUpdate(CView*, LPARAM, CObject*);
	void OnTemperature();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedStart();
	afx_msg void OnBnClickedSequencebutton();
	afx_msg void OnBnClickedSequence2button();
	afx_msg void OnBnClickedSavebutton();
	CString m_seq1;
	CString m_seq2;
	CString m_save;
};


