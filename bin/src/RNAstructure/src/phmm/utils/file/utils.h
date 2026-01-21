#ifndef _UTILS_
#define _UTILS_

#include <stdio.h>
#include "../../../defines.h"

#define __MAX_PATH (10000)

const char* resolve_data_dir();
bool check_file(char* fp);
void validate_file(char* fp);
char* x_fgets(char* buff, int size, FILE* file);
FILE* open_f(const char* fp, const char* mode);

#endif // _UTILS_

