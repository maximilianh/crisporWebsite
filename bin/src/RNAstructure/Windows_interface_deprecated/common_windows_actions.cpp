#include "common_windows_actions.h"


//Get a sequence dialog and return the resulting filename
//Return false if canceled by user or true when it works
//Startpath is a cstring startpath
bool GetSequenceDialog(CString *filename, char *startpath) {

	CFileDialog *filedialog;
	
	filedialog = new CFileDialog(TRUE,NULL,"",OFN_FILEMUSTEXIST|OFN_HIDEREADONLY,
		"Sequence Files (*.seq)|*.seq||");
	
	filedialog->m_ofn.lpstrInitialDir=startpath;
	if (filedialog->DoModal()==IDOK) {
		//strcpy(m_sequencename.GetBuffer(10),(filedialog->GetPathName()).GetBuffer(0));
		(*filename)=(filedialog->GetPathName()).GetBuffer(30);

		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(startpath,_MAX_PATH);

		//now store the path in Startpath so that the program can start here next time:
		//_getcwd(GetFoldDocument()->startpath,_MAX_PATH);
		CString path;
		path = filedialog->GetPathName();
		int i = path.GetLength();
		while(i>=0){
			
			if (path[i]=='\\') break;
			i--;
		}
		if (i>_MAX_PATH) i = _MAX_PATH;
		strncpy(startpath,path.GetBuffer(1),i);
		*(startpath + i) ='\0';

		return true;
	}

	//IDOK wasn't returned by the file dialog; action must have been canceled
	return false;

}

//A function for getting CT names:
//Return false if canceled by user or true when it works
//Startpath is a cstring startpath
bool GetCTDialog(CString *filename, char *startpath) {

	//The user is specifying the CT file name explicitly
	CFileDialog *filedialog;
	filedialog = new CFileDialog(FALSE,".ct",*filename,OFN_OVERWRITEPROMPT|OFN_HIDEREADONLY,
		"CT Files (*.ct)|*.ct|All Files|*.*||");
	if (filedialog->DoModal()==IDOK) {

		*filename=(filedialog->GetPathName()).GetBuffer(0);
			

		delete filedialog;
		return true;

	}
	delete filedialog;
	return false;

}

//Manufacture a CT filename from a .seq filename
void GetCTName(CString *Sequence, CString *CT) {
	char *ctname;
	short int i;
	i = Sequence->GetLength();
	
	ctname = new char[i+4];//allocate enough space so that 
													//three characters can be added 
													//to the name if necessary
	strcpy(ctname,Sequence->GetBuffer(10));
	//count the characters to the .
	
	while(i>=0){
		
		if (ctname[i]=='.') break;
		i--;
	}
	if (i==0) i = Sequence->GetLength();
	strcpy(ctname+i+1,"ct\0");
	*CT=ctname;
	
	delete[] ctname;
		

}