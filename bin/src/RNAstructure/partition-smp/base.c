#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include "base.h"
#include "util.h"

base_t base_from_char(const char c)
{
  switch (c) {
  case 'X':
    return X;
  case 'A':
  case 'a':
    return A;
  case 'C':
  case 'c':
    return C;
  case 'G':
  case 'g':
    return G;
  case 'U':
  case 'u':
  case 'T':
  case 't':
    return U;
  default:
    die("base_from_char: unknown base: \"%c\"", c);
    return A; /* never get here */
  }
}

char base_as_char(const base_t b)
{
  switch (b) {
  case X:
    return 'X';
  case A:
    return 'A';
  case C:
    return 'C';
  case G:
    return 'G';
  case U:
    return 'U';
  default:
    die("base_as_char: unknown base");
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

static int is_valid_base(char c)
{
  return strchr("XAaCcGgUuTt", c) != 0;
}

static int contains_only_valid_bases(const char *s)
{
  while (*s) {
    if (!is_valid_base(*s))
      return 0;
    s++;
  }
  return 1;
}

char *sequence(const char *arg)
{
  return contains_only_valid_bases(arg) ?
    strdup(arg) : sequence_from_file(arg);
}

char *sequence_from_file(const char *fn)
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

void sequence_from_string(base_t *b, const char *s)
{
  for (; *s; s++, b++) {
    *b = base_from_char(*s);
  }
}

#if 0
base_t *sequence_from_string(const char *s)
{
  base_t *b = (base_t *) safe_malloc(strlen(s)*sizeof(base_t));
  int n = 0;
  while (*s)
    b[n++] = base_from_char(*s++);
  return b;
}
#endif
