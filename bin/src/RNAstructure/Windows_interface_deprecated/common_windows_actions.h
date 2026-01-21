#pragma once
#include "stdafx.h"
#include "afxwin.h"


//A function for getting sequence names:
bool GetSequenceDialog(CString *filename, char *startpath);

//A function for getting CT names:
bool GetCTDialog(CString *filename, char *startpath);

//Manufacture a CT filename from a .seq filename
void GetCTName(CString *Sequence, CString *CT);