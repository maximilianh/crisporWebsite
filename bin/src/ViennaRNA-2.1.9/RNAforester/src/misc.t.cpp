#ifndef _MISC_T_CPP_
#define _MISC_T_CPP_

// #include <stdio.h>

#include "misc.h"

template <class T>
void showArray(T *array, Uint m,Uint n)
{
//  char buf[20];
  Uint i,maxElemLen;
  T maxValue;

  // check for the largest value
  maxValue=0;
  for(i=0;i<m*n;i++)
    maxValue = max(maxValue,array[i]);

  //snprintf(buf,20,"%.02f",maxValue);
  //maxElemLen=strlen(buf)+1;
  maxElemLen=5;	

  for(i=0;i<m*n;i++)
    {
    if(i%n==0)
      cout << endl;

    //sprintf(buf,"%*.02f ",maxElemLen,array[i]);
	cout << "[" << i << "]" << array[i] << " ";
//    cout << buf;
    }
  cout << endl;
}


#endif
