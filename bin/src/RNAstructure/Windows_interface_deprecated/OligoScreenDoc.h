#pragma once


// COligoScreenDoc document

class COligoScreenDoc : public CDocument
{
	DECLARE_DYNCREATE(COligoScreenDoc)

public:
	COligoScreenDoc();
	COligoScreenDoc(char *DataPath, char *StartPath);
	CFrameWnd *Frame;
	virtual ~COligoScreenDoc();
	virtual void Serialize(CArchive& ar);   // overridden for document i/o
	char *datapath,*startpath;
	
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:
	virtual BOOL OnNewDocument();

	DECLARE_MESSAGE_MAP()
};
