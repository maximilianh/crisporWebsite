// RNAstructureView.cpp : implementation of the CRNAstructureView class
//

#include "stdafx.h"
#include "RNAstructure.h"

#include "RNAstructureDoc.h"
#include "RNAstructureView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CRNAstructureView

IMPLEMENT_DYNCREATE(CRNAstructureView, CFormView)

BEGIN_MESSAGE_MAP(CRNAstructureView, CFormView)
	//{{AFX_MSG_MAP(CRNAstructureView)
		// NOTE - the ClassWizard will add and remove mapping macros here.
		//    DO NOT EDIT what you see in these blocks of generated code!
	//}}AFX_MSG_MAP
	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, CFormView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, CFormView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, CFormView::OnFilePrintPreview)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CRNAstructureView construction/destruction

CRNAstructureView::CRNAstructureView()
	: CFormView(CRNAstructureView::IDD)
{
	//{{AFX_DATA_INIT(CRNAstructureView)
		// NOTE: the ClassWizard will add member initialization here
	//}}AFX_DATA_INIT
	// TODO: add construction code here

}

CRNAstructureView::~CRNAstructureView()
{
}

void CRNAstructureView::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CRNAstructureView)
		// NOTE: the ClassWizard will add DDX and DDV calls here
	//}}AFX_DATA_MAP
}

BOOL CRNAstructureView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CFormView::PreCreateWindow(cs);
}

void CRNAstructureView::OnInitialUpdate()
{
	CFormView::OnInitialUpdate();
	ResizeParentToFit();

}

/////////////////////////////////////////////////////////////////////////////
// CRNAstructureView printing

BOOL CRNAstructureView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// default preparation
	return DoPreparePrinting(pInfo);
}

void CRNAstructureView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add extra initialization before printing
}

void CRNAstructureView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: add cleanup after printing
}

void CRNAstructureView::OnPrint(CDC* pDC, CPrintInfo* /*pInfo*/)
{
	// TODO: add customized printing code here
}

/////////////////////////////////////////////////////////////////////////////
// CRNAstructureView diagnostics

#ifdef _DEBUG
void CRNAstructureView::AssertValid() const
{
	CFormView::AssertValid();
}

void CRNAstructureView::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}

CRNAstructureDoc* CRNAstructureView::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CRNAstructureDoc)));
	return (CRNAstructureDoc*)m_pDocument;
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CRNAstructureView message handlers
