#ifndef _DEBUG_H_
#define _DEBUG_H_

#include <iostream>

#ifndef RNAF_DEBUGLEVEL
   #define RNAF_DEBUGLEVEL 5 
#endif

#define DBG_OFF 6
#define DBG_QWATCH                1

// program intern debug levels
#define DBG_BACKTRACK             DBG_OFF
#define DBG_GET_PROFILE_STRUCTURE DBG_OFF
#define DBG_ALGEBRA               DBG_OFF
#define DBG_ALIGNMENT             DBG_OFF
#define DBG_MULTIPLE              DBG_OFF

#ifdef NDEBUG
  #define TRACE(L,C,M)
  #define WATCH(L,C,M)
  #define QWATCH(L,C,M)
#else
  #define TRACE(L,C,M)     if(L <= RNAF_DEBUGLEVEL) cout << C << " - " << M << endl;
  #define WATCH(L,C,M)     if(L <= RNAF_DEBUGLEVEL) cout << C << " - " << #M  << ": " << M << endl;
  #define QWATCH(M)        if(DBG_QWATCH <= RNAF_DEBUGLEVEL) cout << #M  << ": " << M << endl;
#endif

#endif
