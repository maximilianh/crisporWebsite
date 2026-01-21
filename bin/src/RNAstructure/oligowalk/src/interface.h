
#if !defined (__INTERFACE)
#define __INTERFACE
#include <iostream>
using namespace std;



///////////////////////////////////////////////////////////////////////////
////////TProgressDialog is opened to notify the user of the progress of
////////the folding algorithm

class TProgressDialog 
{
	public:

   	void update(int frac) {
		cout << frac <<"%"<<flush;
	}


};


/////////////////////////////////////////////////////////////////////////
/////////TRnaDialog is the main program interface


#endif
