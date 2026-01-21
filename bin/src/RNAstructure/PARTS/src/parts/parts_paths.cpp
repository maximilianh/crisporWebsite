#include <stdio.h>
#include "parts_paths.h"
#include "parts_compilation_directives.h"

// Concatenate file path and file name using appropriate delimiter.
void concat_path_fn(char* path_fn, char* path, char* fn)
{
#ifdef WINDOWS_PATHS
	sprintf(path_fn, "%s\\%s", path, fn);
#endif

#ifdef LINUX_PATHS
	sprintf(path_fn, "%s/%s", path, fn);
#endif
}
