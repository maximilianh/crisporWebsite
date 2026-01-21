#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

void die(const char *s, ...);
void *safe_malloc(size_t size);
FILE *safe_fopen(const char *path, const char *mode);
int end_of_file(FILE *f);
int isdigits(const char *);
int string_begins_with(const char *buf, const char *start);
  
#ifdef __cplusplus
}
#endif

#endif /* UTIL_H */
