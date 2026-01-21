#pragma once

#include "oligowalk.h"

// COligoWalkRNAiView form view

class COligoWalkRNAiView : public COligoWalk
{
	DECLARE_DYNCREATE(COligoWalkRNAiView)

protected:
	COligoWalkRNAiView();           // protected constructor used by dynamic creation
	virtual ~COligoWalkRNAiView();

public:
	enum { IDDI = IDD_OLIGOWALK_RNAI };
	void OnInitialUpdate();
	
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
};


