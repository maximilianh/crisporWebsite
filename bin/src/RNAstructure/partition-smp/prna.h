#ifndef PRNA_H
#define PRNA_H

#include "param.h"

#ifdef __cplusplus
extern "C" {
#endif
struct prna
{
  int n; /* number of bases */
  base_t *seq;  /* sequence */
  int *base_can_pair; /* ?? (or is it constrained not to pair) */
  real_t *v; /* n x n array */
  real_t *w5, *w3; /* n elements */
};
  typedef struct prna *prna_t;
  prna_t prna_new(const char *s, param_t par, int quiet, int *base_cp);
  void prna_delete(prna_t);
  void prna_show(const prna_t);
  void prna_write_neg_log10_probabilities(const prna_t, const char *fn);
  void prna_write_probability_matrix(const prna_t, const char *fn);
  void prna_write_probknot(const prna_t, const char *fn, const char *s, int min_helix_length);
  real_t* array_val(real_t *a, int i, int j, int n, const int *bcp);
  real_t probability_of_pair(const prna_t p, int i, int j);
  
  int *generate_bcp(const char *s);

  real_t get_v_array(const prna_t p, int i, int j);
  real_t get_w3_array(const prna_t p, int i);
  real_t get_w5_array(const prna_t p, int i);

#ifdef __cplusplus
}
#endif

#endif /* PRNA_H */
