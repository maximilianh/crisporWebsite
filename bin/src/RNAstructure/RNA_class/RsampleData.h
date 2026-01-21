#ifndef R_SAMPLE_DATA_H
#define R_SAMPLE_DATA_H

#include <vector>
#include "../src/defines.h"

#include <stdlib.h>

using namespace std;

class RsampleData {

  public:
        RsampleData(bool isDMS, double max, const char* const uu_file = NULL, const char* const pend_file = NULL, const char* const  pmid_file = NULL);
        static const char* GetErrorMessage(const int code);
		
	// reactivities for PUE
	vector<double> react_pend;
	vector<double> react_pmid;
	vector<double> react_uu;
	
	int ErrorCode;
};
#endif // R_SAMPLE_DATA_H
