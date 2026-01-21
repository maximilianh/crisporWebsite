#pragma once


// CColorKeyDoc document

class CColorKeyDoc : public CDocument
{
	DECLARE_DYNCREATE(CColorKeyDoc)

public:
	CColorKeyDoc();
	virtual ~CColorKeyDoc();
	virtual void Serialize(CArchive& ar);   // overridden for document i/o
	
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:
	virtual BOOL OnNewDocument();

	DECLARE_MESSAGE_MAP()
};
