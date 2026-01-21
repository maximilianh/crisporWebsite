#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include "fbase.h"
#include "../partition-smp/util.h"

fbase_t fbase_from_char(const char c)
{
  switch (c) {
  case 'A':
    return A;
  case 'C':
    return C;
  case 'G':
    return G;
  case 'U':
  case 'T':
    return U;
  default:
    die("fbase_from_char: unknown base: \"%c\"", c);
    return A; /* never get here */
  }
}

char fbase_as_char(const fbase_t b)
{
  switch (b) {
  case A:
    return 'A';
  case C:
    return 'C';
  case G:
    return 'G';
  case U:
    return 'U';
  default:
    die("fbase_as_char: unknown base");
    return 0;
  }
}

static void read_line(FILE *f)
{
  char c;
  do { c = fgetc(f); } while (c != '\n' && c != EOF);
}

static void skip_header(FILE *f)
{
  int extra = 0;
  char c;
  for (c = fgetc(f); c == ';' || c == '>'; read_line(f), c = fgetc(f))
    if (c == ';')
      extra = 1;
    else
      extra = 0;
  ungetc(c,f);
  if (extra)
    read_line(f);
}

static int contains_only_valid_bases(const char *s)
{
  while (*s == 'A' || *s == 'C' || *s == 'G' || *s == 'U' || *s == 'T')
    s++;
  return !*s;
}

char *fsequence(const char *arg)
{
  return contains_only_valid_bases(arg) ?
    strdup(arg) : fsequence_from_file(arg);
}

char *fsequence_from_file(const char *fn)
{
  char *s = 0;
  int n;
  char c;
  FILE *f = safe_fopen(fn, "r");
  skip_header(f);
  n = 0;
  while ((c = fgetc(f)), c != EOF && c != '1')
    if (!isspace(c))
      n++;
  s = (char *) safe_malloc(n+1);
  rewind(f);
  skip_header(f);
  n = 0;
  while ((c = fgetc(f)), c != EOF && c != '1')
    if (!isspace(c))
      s[n++] = c;
  s[n] = 0;
  fclose(f);
  return s;
}

fbase_t *fsequence_from_string(const char *s)
{
  fbase_t *b = (fbase_t *) safe_malloc(strlen(s)*sizeof(fbase_t));
  int n = 0;
  while (*s)
    b[n++] = fbase_from_char(*s++);
  return b;
}

