//{{AFX_DATA_MAP(COligoWalk)
	DDX_Text(pDX, IDC_LENGTH, m_length);
	DDV_MinMaxInt(pDX, m_length, 1, 100);
	DDX_Text(pDX, IDC_START, m_start);
	DDX_Text(pDX, IDC_STOP, m_stop);
	DDX_Text(pDX, IDC_CTNAME, m_ctname);
	DDX_Text(pDX, IDC_REPORTNAME, m_reportname);
	DDX_Text(pDX, IDC_CONCENTRATION, m_concentration);
	DDV_MinMaxDouble(pDX, m_concentration, 0., 100000.);
	//}}AFX_DATA_MAP