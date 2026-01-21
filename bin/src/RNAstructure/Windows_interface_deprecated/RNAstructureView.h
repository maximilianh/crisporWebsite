// RNAstructureView.h : interface of the CRNAstructureView class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_RNASTRUCTUREVIEW_H__56CD85EE_8D8A_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_RNASTRUCTUREVIEW_H__56CD85EE_8D8A_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


class CRNAstructureView : public CFormView
{
protected: // create from serialization only
	CRNAstructureView();
	DECLARE_DYNCREATE(CRNAstructureView)

public:
	//{{AFX_DATA(CRNAstructureView)
	enum{ IDD = IDD_RNASTRUCTURE_FORM };
		// NOTE: the ClassWizard will add data members here
	//}}AFX_DATA

// Attributes
public:
	CRNAstructureDoc* GetDocument();

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CRNAstructureView)
	public:
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	virtual void OnInitialUpdate(); // called first time after construct
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnPrint(CDC* pDC, CPrintInfo* pInfo);
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CRNAstructureView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	//{{AFX_MSG(CRNAstructureView)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in RNAstructureView.cpp
inline CRNAstructureDoc* CRNAstructureView::GetDocument()
   { return (CRNAstructureDoc*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_RNASTRUCTUREVIEW_H__56CD85EE_8D8A_11D4_9F32_00C0F02A5F5D__INCLUDED_)
