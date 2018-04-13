#ifndef _MISC_H
#define _MISC_H

#include "types.h"

#define DELETE(T)   if(T) \
                    { \
                      delete T; \
                      T=NULL; \
                    }

#define DELETE_ARRAY(T)  if(T) \
                         { \
                           delete[] T; \
                           T=NULL; \
                         }

template <class T>
void showArray(T *array, Uint m,Uint n);

#endif
