#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdarg.h>
#include <ctype.h>
#include "util.h"

void die(const char *fmt, ...)
{
  va_list argp;
  va_start(argp, fmt);
  char *ret;
  vasprintf(&ret, fmt, argp);
  va_end(argp);
  fprintf(stderr, "%s\n", ret);
  free(ret);
  exit(1);
}

void *safe_malloc(size_t size)
{
  void *p = malloc(size);
  if (!p)
    die("safe_malloc: could not allocate %lu bytes", size);
  return p;
}

FILE *safe_fopen(const char *path, const char *mode)
{
  FILE *f = fopen(path,mode);
  if (!f)
    die("safe_fopen: could not open file \'%s\'", path);
  return f;
}

int end_of_file(FILE *f)
{
  int c = fgetc(f);
  if (c == EOF)
    return 1;
  ungetc(c, f);
  return 0;
}

int isdigits(const char *s)
{
  for ( ; *s; s++)
    if (!isdigit(*s))
      return 0;
  return 1;
}

int string_begins_with(const char *buf, const char *start)
{
  return !strncmp(buf, start, strlen(start));
}
