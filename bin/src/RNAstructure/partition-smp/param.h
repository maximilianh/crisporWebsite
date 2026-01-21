#ifndef PARAM_H
#define PARAM_H

#include "base.h"
#include "real.h"

#ifdef __cplusplus
extern "C" {
#endif

extern const double T37;// = RCONST(310.15);
extern const double R;// = RCONST(1.987213e-3);//gas constant
extern const double RT;// = RCONST(1.987213e-3)*RCONST(310.15);;
#define LOOP_MIN 3
#define LOOP_MAX 30
#define LOOKUP_TABLE_MAX 120

typedef real_t tab3_t[NBASE][NBASE][NBASE];
typedef real_t tab4_t[NBASE][NBASE][NBASE][NBASE];
typedef real_t tab6_t[NBASE][NBASE][NBASE][NBASE][NBASE][NBASE];
typedef real_t tab7_t[NBASE][NBASE][NBASE][NBASE][NBASE][NBASE][NBASE];
//typedef real_t tab8_t[NBASE][NBASE][NBASE][NBASE][NBASE][NBASE][NBASE][NBASE];
typedef real_t tab8_t[5][5][5][5][5][5][5][5];

typedef struct param
{
  int use_dna_params;
  int use_enthalpy_params;
  real_t temperature;

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

  real_t internal_loop_initiation[LOOP_MAX+1];
  real_t bulge_loop_initiation[LOOP_MAX+1];
  real_t hairpin_loop_initiation[LOOP_MAX+1];
  real_t Extrapolation_for_large_loops;
  real_t maximum_correction;
  real_t fm_array_first_element;
  real_t a, multibranched_loop_offset;
  real_t b, multibranched_loop_per_nuc_penalty;
  real_t c, multibranched_loop_helix_penalty;
  real_t a_2c, a_2b_2c;
  real_t terminal_AU_penalty;
  real_t bonus_for_GGG_hairpin;
  real_t c_hairpin_slope;
  real_t c_hairpin_intercept;
  real_t c_hairpin_of_3;
  real_t Bonus_for_Single_C_bulges_adjacent_to_C;

  int ntriloop;
  struct triloop_t { base_t seq[5]; real_t val; } triloop[LOOKUP_TABLE_MAX];
  int ntloop;
  struct tloop_t { base_t seq[6]; real_t val; } tloop[LOOKUP_TABLE_MAX];
  int nhexaloop;
  struct hexaloop_t { base_t seq[8]; real_t val; } hexaloop[LOOKUP_TABLE_MAX];

  tab3_t dangle_3p;
  tab3_t dangle_5p;

  tab6_t int11;
  tab8_t int22;
  tab7_t int21;
  
} *param_t;

void param_read_from_text(const char *path, param_t p, int use_dna_params, int use_enthapy_params);
void param_read_from_binary(const char *path, param_t p);
void param_save_to_binary(const char *path, const param_t p);
void param_show(const param_t p);
void param_read_alternative_temperature(const char *path, param_t p, 
                                        real_t temperature, int use_dna_params);
#ifdef __cplusplus
}
#endif

#endif /* PARAM_H */
