#ifndef BASE_H
#define BASE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "cu.h"

#define NBASE 5
typedef enum { A, C, G, U } base_t; 
char base_as_char(const base_t);
base_t base_from_char(const char);
char *sequence(const char *arg);
char *sequence_from_file(const char *fname);
base_t *sequence_from_string(const char *s);

DEV HOST inline static int is_cp(const base_t i, const base_t j)
{
  return (i == A && j == U) || (i == C && j == G) || (i == G && j == U);
}

DEV HOST inline static int is_canonical_pair(const base_t i, const base_t j)
{
  return is_cp(i,j) || is_cp(j,i);
}

DEV inline static int contains_only_base(const base_t b, const int n, const base_t *seq)
{
  int i;
  for (i = 0; i < n; i++)
    if (seq[i] != b)
      return 0;
  return 1;
}

DEV inline static int sequences_match(const base_t *s1, const base_t *s2, int n)
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
