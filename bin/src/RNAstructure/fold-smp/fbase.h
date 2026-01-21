#ifndef BASE_H
#define BASE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "cu.h"

#define NBASE 4
typedef enum { A, C, G, U } fbase_t; 
char fbase_as_char(const fbase_t);
fbase_t fbase_from_char(const char);
char *fsequence(const char *arg);
char *fsequence_from_file(const char *fname);
fbase_t *fsequence_from_string(const char *s);

DEV HOST inline static int is_cp(const fbase_t i, const fbase_t j)
{
  return (i == A && j == U) || (i == C && j == G) || (i == G && j == U);
}

DEV HOST inline static int is_canonical_pair(const fbase_t i, const fbase_t j)
{
  return is_cp(i,j) || is_cp(j,i);
}

DEV inline static int contains_only_base(const fbase_t b, const int n, const fbase_t *seq)
{
  int i;
  for (i = 0; i < n; i++)
    if (seq[i] != b)
      return 0;
  return 1;
}

DEV inline static int sequences_match(const fbase_t *s1, const fbase_t *s2, int n)
{
  int i;
  for (i = 0; i < n; i++)
    if (s1[i] != s2[i])
      return 0;
  return 1;
}
 
#ifdef __cplusplus
}
#endif

#endif /* BASE_H */
