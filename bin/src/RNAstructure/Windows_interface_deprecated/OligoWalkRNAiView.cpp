// OligoWalkRNAiView.cpp : implementation file
//

#include "stdafx.h"
#include "RNAstructure.h"
#include "OligoWalkRNAiView.h"


// COligoWalkRNAiView

IMPLEMENT_DYNCREATE(COligoWalkRNAiView, CFormView)

COligoWalkRNAiView::COligoWalkRNAiView()
	:COligoWalk(COligoWalkRNAiView::IDDI) 
{
	m_asuf = 0;
	m_tofe = -10;
	m_fnnfe = -1.33;

}

COligoWalkRNAiView::~COligoWalkRNAiView()
{
}

void COligoWalkRNAiView::OnInitialUpdate() {

	COligoWalk::OnInitialUpdate();
}
void COligoWalkRNAiView::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
#include "commonoligowalkdataexchange.h"
	DDX_Text(pDX, IDC_ASUF, m_asuf);
	DDX_Text(pDX, IDC_TOFE, m_tofe);
	DDX_Text(pDX, IDC_FNNFE, m_fnnfe);
}

BEGIN_MESSAGE_MAP(COligoWalkRNAiView, CFormView)
#include "commonoligowalkmessagemap.h"	
END_MESSAGE_MAP()


// COligoWalkRNAiView diagnostics

#ifdef _DEBUG
void COligoWalkRNAiView::AssertValid() const
{
	CFormView::AssertValid();
}

void COligoWalkRNAiView::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}
#endif //_DEBUG


// COligoWalkRNAiView message handlers
