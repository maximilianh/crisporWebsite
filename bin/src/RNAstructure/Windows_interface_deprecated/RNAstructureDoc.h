// RNAstructureDoc.h : interface of the CRNAstructureDoc class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_RNASTRUCTUREDOC_H__56CD85EC_8D8A_11D4_9F32_00C0F02A5F5D__INCLUDED_)
#define AFX_RNASTRUCTUREDOC_H__56CD85EC_8D8A_11D4_9F32_00C0F02A5F5D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


class CRNAstructureDoc : public CDocument
{
protected: // create from serialization only
	CRNAstructureDoc();
	DECLARE_DYNCREATE(CRNAstructureDoc)

// Attributes
public:

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CRNAstructureDoc)
	public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CRNAstructureDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	//{{AFX_MSG(CRNAstructureDoc)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_RNASTRUCTUREDOC_H__56CD85EC_8D8A_11D4_9F32_00C0F02A5F5D__INCLUDED_)
