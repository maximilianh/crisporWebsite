/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#include "errmsg.h"

#include <iostream>
#include <cstdlib>

using namespace std;

// Function for outputting info in case of an error
void errmsg(int err,int erri) {
  if (err==30) {
    cout << "End Reached at traceback #"<<erri<<"\n";
    exit(1);
  }
  if (err==100) {
    cout << "error # "<<erri;
    exit(1);
  }
  switch (err) {
	case 1:
   	cout << "Could not allocate enough memory";
    break;
  case 2:
   	cout << "Too many possible base pairs";
    break;
  case 3:
   	cout << "Too many helixes in multibranch loop";
  case 4:
   	cout << "Too many structures in CT file";
  default:
   	cout << "Unknown error";
  }
  cin >> err;
  exit(1);
  return;
}
