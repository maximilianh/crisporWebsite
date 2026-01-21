#ifndef PARAM_H
#define PARAM_H

#include "fbase.h"
#include "int.h"

#ifdef __cplusplus
extern "C" {
#endif

#define RT (RCONST(1.987213e-3)*RCONST(310.15))
#define LOOP_MIN 3
#define LOOP_MAX 30
#define LOOKUP_TABLE_MAX 120


typedef int_t tab3_t[NBASE][NBASE][NBASE];
typedef int_t tab4_t[NBASE][NBASE][NBASE][NBASE];
typedef int_t tab6_t[NBASE][NBASE][NBASE][NBASE][NBASE][NBASE];
typedef int_t tab7_t[NBASE][NBASE][NBASE][NBASE][NBASE][NBASE][NBASE];
typedef int_t tab8_t[NBASE][NBASE][NBASE][NBASE][NBASE][NBASE][NBASE][NBASE];

typedef struct fparam
{
  int use_dna_fparams;

  tab4_t coaxial;
  tab4_t coaxstack;
  tab4_t stack;
  tab4_t tstack;
  tab4_t tstackcoax;
  tab4_t tstackh;
  tab4_t tstacki;
  tab4_t tstacki1n;
  tab4_t tstacki23;
  tab4_t tstackm;

  int_t internal_loop_initiation[LOOP_MAX+1];
  int_t bulge_loop_initiation[LOOP_MAX+1];
  int_t hairpin_loop_initiation[LOOP_MAX+1];
  int_t Extrapolation_for_large_loops;
  float prelog;
  int_t maximum_correction;
  int_t fm_array_first_element;
  int_t a, multibranched_loop_offset;
  int_t b, multibranched_loop_per_nuc_penalty;
  int_t c, multibranched_loop_helix_penalty;
  int_t a_2c, a_2b_2c;
  int_t terminal_AU_penalty;
  int_t bonus_for_GGG_hairpin;
  int_t c_hairpin_slope;
  int_t c_hairpin_intercept;
  int_t c_hairpin_of_3;
  int_t Bonus_for_Single_C_bulges_adjacent_to_C;

  int ntriloop;
  struct triloop_t { fbase_t seq[5]; int_t val; } triloop[LOOKUP_TABLE_MAX];
  int ntloop;
  struct tloop_t { fbase_t seq[6]; int_t val; } tloop[LOOKUP_TABLE_MAX];
  int nhexaloop;
  struct hexaloop_t { fbase_t seq[8]; int_t val; } hexaloop[LOOKUP_TABLE_MAX];

  tab3_t dangle_3p;
  tab3_t dangle_5p;

  tab6_t int11;
  tab8_t int22;
  tab7_t int21;
  
} *fparam_t;

void fparam_read_from_text(const char *path, int use_dna_fparams, fparam_t p);
void fparam_read_from_binary(const char *path, fparam_t p);
void fparam_save_to_binary(const char *path, const fparam_t p);
void fparam_show(const fparam_t p);

#ifdef __cplusplus
}
#endif

#endif /* PARAM_H */
