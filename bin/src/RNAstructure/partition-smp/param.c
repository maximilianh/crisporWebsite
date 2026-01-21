#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <ctype.h>
#include "param.h"
#include "util.h"

const double T37 = RCONST(310.15);
const double R = RCONST(1.987213e-3); /* gas constant */
const double RT = RCONST(1.987213e-3)*RCONST(310.15);

real_t init_val(int *b, int n, real_t Xval)
{
  int i;
  for (i = 0; i < n; i++)
    if (b[i] == X)
      return Xval;
  return NOT_A_NUMBER;
}

static void init_tab3(tab3_t *t, real_t Xval)
{
  int i[3];
  for (i[0] = 0; i[0] < NBASE; i[0]++)
    for (i[1] = 0; i[1] < NBASE; i[1]++)
      for (i[2] = 0; i[2] < NBASE; i[2]++)
	(*t)[i[0]][i[1]][i[2]] = init_val(i,3,Xval);
}

static void init_tab4(tab4_t *t, real_t Xval)
{
  int i[4];
  for (i[0] = 0; i[0] < NBASE; i[0]++)
    for (i[1] = 0; i[1] < NBASE; i[1]++)
      for (i[2] = 0; i[2] < NBASE; i[2]++)
	for (i[3] = 0; i[3] < NBASE; i[3]++)
	  (*t)[i[0]][i[1]][i[2]][i[3]] = init_val(i,4,Xval);
}

static void init_tab6(tab6_t *t, real_t Xval)
{
  int i[6];
  for (i[0] = 0; i[0] < NBASE; i[0]++)
    for (i[1] = 0; i[1] < NBASE; i[1]++)
      for (i[2] = 0; i[2] < NBASE; i[2]++)
	for (i[3] = 0; i[3] < NBASE; i[3]++)
	  for (i[4] = 0; i[4] < NBASE; i[4]++)
	    for (i[5] = 0; i[5] < NBASE; i[5]++)
	      (*t)[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]] = init_val(i,6,Xval);
}

static void init_tab7(tab7_t *t, real_t Xval)
{
  int i[7];
  for (i[0] = 0; i[0] < NBASE; i[0]++)
    for (i[1] = 0; i[1] < NBASE; i[1]++)
      for (i[2] = 0; i[2] < NBASE; i[2]++)
	for (i[3] = 0; i[3] < NBASE; i[3]++)
	  for (i[4] = 0; i[4] < NBASE; i[4]++)
	    for (i[5] = 0; i[5] < NBASE; i[5]++)
	      for (i[6] = 0; i[6] < NBASE; i[6]++)
		(*t)[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]][i[6]] = init_val(i,7,Xval);
}

static void init_tab8(tab8_t *t, real_t Xval)
{
  int i[8];
  for (i[0] = 0; i[0] < NBASE; i[0]++)
    for (i[1] = 0; i[1] < NBASE; i[1]++)
      for (i[2] = 0; i[2] < NBASE; i[2]++)
	for (i[3] = 0; i[3] < NBASE; i[3]++)
	  for (i[4] = 0; i[4] < NBASE; i[4]++)
	    for (i[5] = 0; i[5] < NBASE; i[5]++)
	      for (i[6] = 0; i[6] < NBASE; i[6]++)
		for (i[7] = 0; i[7] < NBASE; i[7]++)
		  (*t)[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]][i[6]][i[7]] = init_val(i,8,Xval);
}

#define MAXLINE 1024

FILE *parfile(const char *name, int use_dna_params, int use_enthalpy_params)
{
  char buf[MAXLINE+1];
  const char *prefix = use_dna_params ? "dna" : "rna";
  const char *extension = use_enthalpy_params ? "dh" : "dg";
  if (strlen(name) + 8 > MAXLINE)
    die("parfile: name too long: %s", name);
  sprintf(buf, "%s.%s.%s", prefix, name, extension);
  return safe_fopen(buf,"r");
}

static void look_for_line_containing(FILE *f, const char *s)
{
  char buf[MAXLINE+1];
  while (fgets(buf, MAXLINE, f))
    if (strstr(buf, s))
      return;
  die("look_for_line_starting_with: couldn't find %s", s);
}

static const char whitespace[] = " \f\n\r\t\v";

static int is_only_whitespace(const char *buf)
{
  const char *b;
  for (b = buf; *b; b++)
    if (!isspace(*b))
      return 0;
  return 1;
}

static int read_next_line(FILE *f, char *buf)
{
  if (!fgets(buf, MAXLINE, f))
    return 0;
  if (strlen(buf) >= MAXLINE)
    die("read_next_line: line too long");
  if (buf[0] == '#')
    return read_next_line(f, buf); /* comment */
  if (is_only_whitespace(buf))
    return read_next_line(f, buf); /* empty line */
  return 1;
}

static void expect_line_containing(FILE *f, const char *s)
{
  char buf[MAXLINE+1];
  if (!(read_next_line(f, buf) && strstr(buf, s)))
    die("expect_line_containing: did not find %s", s);
}

static void look_for_dashes(FILE *f)
{
  look_for_line_containing(f, "--------");
}

static void look_for_arrow(FILE *f)
{
  look_for_line_containing(f, "-->");
}

static real_t value_from_string(const char *s)
{
  return strcmp(s,".") ? STR_TO_REAL(s,0) / RT : NOT_A_NUMBER;
}

#define sqNBASEm1 ((NBASE-1)*(NBASE-1))
typedef char field_t[sqNBASEm1+1][MAXLINE+1];

static int split(char *buf, field_t *f)
{
  char *b;
  int n = 0;
  for (b = strtok(buf, whitespace); b; b = strtok(0, whitespace))
    strcpy(&(*f)[n++][0], b);
  return n;
}

static void split_expecting(char *buf, field_t *f, int nf_expected)
{
  const int nf = split(buf, f);
  if (nf != nf_expected)
    die("split_expecting: found %d fields, expected %d", nf, nf_expected);
}

static void read_stack(const char *name, int use_dna_params, tab4_t *p, int use_enthalpy_params)
{
  init_tab4(p,0);
  FILE *f = parfile(name, use_dna_params, use_enthalpy_params);
  const int coaxial = !strcmp(name,"coaxial");
  char buf[MAXLINE+1];
  field_t fld;
  base_t i, j, y[NBASE];
  while (read_next_line(f, buf)) {
    const int nf = split(buf, &fld);
    if (nf == 1) {
      if (fld[0][1] == 'X')
	i = base_from_char(fld[0][0]);
      else if (fld[0][1] == 'Y')
	j = base_from_char(fld[0][0]);
      else
	die("read_stack: expecting BX or BY, found %s", buf);
    } else if (nf == NBASE-1) {
      int ib;
      for (ib = 0; ib < NBASE-1; ib++)
	y[ib] = base_from_char(fld[ib][0]);
    } else if (nf == NBASE) {
      base_t k = base_from_char(fld[0][0]);
      int ib;
      for (ib = 0; ib < NBASE-1; ib++) {
	base_t l = y[ib];
	real_t val = value_from_string(fld[ib+1]);
	if (coaxial)
	  (*p)[j][i][k][l] = val;
	else
	  (*p)[i][j][k][l] = val;
      }
    } else {
      die("read_stack: found %d fields (expecting 1, %d, or %d) in %s", nf, NBASE-1, NBASE, buf);
    }
  }
  fclose(f);
}

static void read_loop(param_t p)
{
  FILE *f = parfile("loop", p->use_dna_params, p->use_enthalpy_params);
  look_for_dashes(f);
  char buf[MAXLINE+1];
  int i;
  for (i = 0; i <= LOOP_MAX; i++) {
    p->internal_loop_initiation[i] = NOT_A_NUMBER;
    p->bulge_loop_initiation[i] = NOT_A_NUMBER;
    p->hairpin_loop_initiation[i] = NOT_A_NUMBER;
  }
  while (fgets(buf, MAXLINE, f)) {
    if (strlen(buf) >= MAXLINE)
      die("read_loop: line too long");
    char fld[4][MAXLINE];
    if (sscanf(buf, "%s%s%s%s", fld[0], fld[1], fld[2], fld[3]) != 4)
      die("read_loop: could not find 4 fields");
    const int i = atoi(fld[0]);
    if (i <= 0 || i > LOOP_MAX)
      die("read_loop: index out of range");
    p->internal_loop_initiation[i] = value_from_string(fld[1]);
    p->bulge_loop_initiation[i] = value_from_string(fld[2]);
    p->hairpin_loop_initiation[i] = value_from_string(fld[3]);
  }
  fclose(f);
}

static real_t next_value(FILE *f)
{
  look_for_arrow(f);
  char buf[MAXLINE+1];
  read_next_line(f, buf);
  return value_from_string(strtok(buf,whitespace));
}

static void read_next_values(FILE *f, real_t *v, int n)
{
  char buf[MAXLINE+1];
  read_next_line(f, buf);
  int i;
  for (i = 0; i < n; i++)
    v[i] = value_from_string(strtok(i == 0 ? buf : 0, whitespace));
}

static void read_miscloop(param_t p)
{
  FILE *f = parfile("miscloop", p->use_dna_params,p->use_enthalpy_params);
  p->Extrapolation_for_large_loops = next_value(f);
  p->maximum_correction = next_value(f);
  p->fm_array_first_element = next_value(f);
  real_t tmp[3];
  look_for_arrow(f);
  read_next_values(f, tmp, 3);
  p->a = p->multibranched_loop_offset = tmp[0];
  p->b = p->multibranched_loop_per_nuc_penalty = tmp[1];
  p->c = p->multibranched_loop_helix_penalty = tmp[2];
  p->a_2c = p->a + 2*p->c;
  p->a_2b_2c = p->a + 2*p->b + 2*p->c;
  look_for_arrow(f); /* skip efn2 params */
  look_for_arrow(f); /* skip multiloop asym */
  look_for_arrow(f); /* skip multiloop strain */
  p->terminal_AU_penalty = next_value(f);
  p->bonus_for_GGG_hairpin = next_value(f);
  p->c_hairpin_slope = next_value(f);
  p->c_hairpin_intercept = next_value(f);
  p->c_hairpin_of_3 = next_value(f);
  look_for_arrow(f); /* skip Intermolecular initiation */
  p->Bonus_for_Single_C_bulges_adjacent_to_C = next_value(f);
  fclose(f);
}

static void read_small_loop(char *buf, base_t seq[], real_t *val, int n)
{
  char *b = strtok(buf, whitespace);
  int i;
  for (i = 0; i < n; i++)
    seq[i] = base_from_char(b[i]);
  b = strtok(0, whitespace);
  *val = value_from_string(b);
}

#define READ_SMALL_LOOP(X,N)						\
  static void read_##X##loop(param_t p)					\
  {									\
    int n = 0;								\
    FILE *f = parfile(#X"loop", p->use_dna_params,p->use_enthalpy_params);			\
    look_for_dashes(f);							\
    char buf[MAXLINE+1];						\
    while (fgets(buf, MAXLINE, f)) {					\
      if (strlen(buf) >= MAXLINE)					\
	die("read_triloop: line too long");				\
      if (is_only_whitespace(buf))					\
	continue;							\
      read_small_loop(buf, p->X##loop[n].seq, &p->X##loop[n].val, N);	\
      n++;								\
      if (n > LOOKUP_TABLE_MAX)						\
	die("read_"#X"loop: more than %d entries, need to increase LOOKUP_TABLE_MAX", LOOKUP_TABLE_MAX); \
    }									\
    p->n##X##loop = n;							\
    fclose(f);								\
  }

READ_SMALL_LOOP(tri,5)
READ_SMALL_LOOP(t,6)
READ_SMALL_LOOP(hexa,8)

#undef READ_SMALL_LOOP

static void read_dangle(int use_dna_params, int use_enthalpy_params, tab3_t *d3p, tab3_t *d5p)
{
  init_tab3(d3p,0);
  init_tab3(d5p,0);
  FILE *f = parfile("dangle", use_dna_params, use_enthalpy_params);
  char buf[MAXLINE+1];
  base_t i, j, x[NBASE];
  tab3_t *d = 0;
  while (read_next_line(f,buf)) {
    field_t fld0;
    const int nf = split(buf, &fld0);
    if (nf == 1) {
      field_t fld1;
      read_next_line(f,buf);
      split_expecting(buf, &fld1, 1);
      i = base_from_char(fld0[0][0]);
      j = base_from_char(fld1[0][0]);
      if (fld0[0][1] == 'X') {
	d = d3p;
      } else {
	if (fld1[0][1] != 'X')
	  die("read_dangle: expecting BX, found something else", buf);
	d = d5p;
      }
    } else if (nf == NBASE-1) {
      int ib;
      if (strchr(buf, '.'))
	for (ib = 0; ib < NBASE-1; ib++) {
	  base_t k = x[ib];
	  real_t val = value_from_string(fld0[ib]);
	  (*d)[i][j][k] = val;
	}
      else
	for (ib = 0; ib < NBASE-1; ib++) {
	  x[ib] = base_from_char(fld0[ib][0]);
	}
    }
  }
  fclose(f);
}

static void read_int11(int use_dna_params, int use_enthalpy_params, tab6_t *t)
{
  init_tab6(t,INF);
  FILE *f = parfile("int11", use_dna_params, use_enthalpy_params);
  char buf[MAXLINE+1];
  field_t fld;
  while (read_next_line(f, buf)) {
    if (!strstr(buf,"X"))
      die("read_int11: expecting 'X'");
    read_next_line(f, buf);
    split_expecting(buf, &fld, 2);
    const base_t i = base_from_char(fld[0][0]);
    const base_t k = base_from_char(fld[1][0]);
    read_next_line(f, buf);
    split_expecting(buf, &fld, 2);
    const base_t l = base_from_char(fld[1][0]);
    const base_t n = base_from_char(fld[0][0]);
    expect_line_containing(f, "Y");
    read_next_line(f, buf);
    if (split(buf, &fld) != NBASE-1)
      die("read_int11: expecting %d bases", NBASE-1);
    base_t y[NBASE];
    int ib;
    for (ib = 0; ib < NBASE-1; ib++)
      y[ib] = base_from_char(fld[ib][0]);
    for (ib = 0; ib < NBASE-1; ib++) {
      read_next_line(f, buf);
      split_expecting(buf, &fld, NBASE);
      const base_t j = base_from_char(fld[0][0]);
      int jb;
      for (jb = 0; jb < NBASE-1; jb++) {
	const base_t m = y[jb];
	const real_t val = value_from_string(fld[jb+1]);
	(*t)[i][j][k][l][m][n] = val;
      }
    }
  }
  fclose(f);
}

static void read_int22(int use_dna_params, int use_enthalpy_params, tab8_t *t)
{
//  printf("Size : %lu\n",sizeof(*t));
//  printf("Size : %lu\n",sizeof(*t)/sizeof((*t)[0]));
//  printf("Size : %lu\n",sizeof(*t)/sizeof((*t)[0][0]));
//  printf("Size : %lu\n",sizeof(*t)/sizeof((*t)[0][0][0]));
//  printf("Size : %lu\n",sizeof(*t)/sizeof((*t)[0][0][0][0]));
//  printf("Size : %lu\n",sizeof(*t)/sizeof((*t)[0][0][0][0][0]));
//  printf("Size : %lu\n",sizeof(*t)/sizeof((*t)[0][0][0][0][0][0]));
//  printf("Size : %lu\n",sizeof(*t)/sizeof((*t)[0][0][0][0][0][0][0]));
//  printf("Size : %lu\n",sizeof(*t)/sizeof((*t)[0][0][0][0][0][0][0][0]));

//  t = (tab8_t *) safe_malloc(sizeof(tab8_t));

/*
  (*t)[0][0][0][0][0][0][0][0]=0;
  (*t)[1][0][0][0][0][0][0][0]=0;
  (*t)[0][1][0][0][0][0][0][0]=0;
  (*t)[0][0][1][0][0][0][0][0]=0;
  (*t)[0][0][0][1][0][0][0][0]=0;
  (*t)[0][0][0][0][1][0][0][0]=0;
  (*t)[0][0][0][0][0][1][0][0]=0;
  (*t)[0][0][0][0][0][0][1][0]=0;
  (*t)[0][0][0][0][0][0][0][1]=0;

  (*t)[2][0][0][0][0][0][0][0]=0;
  (*t)[0][2][0][0][0][0][0][0]=0;
  (*t)[0][0][2][0][0][0][0][0]=0;
  (*t)[0][0][0][2][0][0][0][0]=0;
  (*t)[0][0][0][0][2][0][0][0]=0;
  (*t)[0][0][0][0][0][2][0][0]=0;
  (*t)[0][0][0][0][0][0][2][0]=0;
  (*t)[0][0][0][0][0][0][0][2]=0;

  (*t)[0][3][0][0][0][0][0][0]=0;
  (*t)[0][0][3][0][0][0][0][0]=0;
  (*t)[0][0][0][3][0][0][0][0]=0;
  (*t)[0][0][0][0][3][0][0][0]=0;
  (*t)[0][0][0][0][0][3][0][0]=0;
  (*t)[0][0][0][0][0][0][3][0]=0;
  (*t)[0][0][0][0][0][0][0][3]=0;
  (*t)[2][1][0][0][0][0][0][0]=0;
  (*t)[2][2][0][0][0][0][0][0]=0;
  (*t)[2][3][0][0][0][0][0][0]=0;
  (*t)[2][4][0][0][0][0][0][0]=0;
  (*t)[2][4][1][0][0][0][0][0]=0;
  (*t)[2][4][1][1][0][0][0][0]=0;
  (*t)[2][4][1][2][0][0][0][0]=0;
  (*t)[2][4][1][3][0][0][0][0]=0;
  (*t)[2][4][1][4][0][0][0][0]=0;
  (*t)[2][4][2][0][0][0][0][0]=0;
  (*t)[2][4][3][0][0][0][0][0]=0;
  (*t)[2][4][4][0][0][0][0][0]=0;
  (*t)[2][4][4][1][0][0][0][0]=0;
  (*t)[2][4][4][2][0][0][0][0]=0;
  (*t)[2][4][4][3][0][0][0][0]=0;
  (*t)[2][4][4][4][0][0][0][0]=0;
  (*t)[2][4][4][4][1][0][0][0]=0;
  (*t)[2][4][4][4][2][0][0][0]=0;
  (*t)[2][4][4][4][3][0][0][0]=0;
  (*t)[2][4][4][4][4][0][0][0]=0;
  (*t)[3][0][0][0][0][0][0][0]=0;

  (*t)[4][4][4][4][4][4][4][4]=0;
*/

  init_tab8(t,INF);
  FILE *f = parfile("int22", use_dna_params, use_enthalpy_params);
  char buf[MAXLINE+1];
  field_t fld;
  base_t y1[sqNBASEm1], y2[sqNBASEm1];
  while (read_next_line(f, buf)) {
    split_expecting(buf, &fld, 4);
    if (strcmp(fld[1], "X1"))
      die("read_int22: expecting 'X1'");
    if (strcmp(fld[2], "Y1"))
      die("read_int22: expecting 'Y1'");
    const base_t i = base_from_char(fld[0][0]);
    const base_t j = base_from_char(fld[3][0]);
    read_next_line(f, buf);
    split_expecting(buf, &fld, 4);
    if (strcmp(fld[1], "X2"))
      die("read_int22: expecting 'X2'");
    if (strcmp(fld[2], "Y2"))
      die("read_int22: expecting 'Y2'");
    const base_t k = base_from_char(fld[0][0]);
    const base_t l = base_from_char(fld[3][0]);
    read_next_line(f, buf);
    split_expecting(buf, &fld, sqNBASEm1);
    int ib;
    for (ib = 0; ib < sqNBASEm1; ib++)
      y1[ib] = base_from_char(fld[ib][0]);
    read_next_line(f, buf);
    split_expecting(buf, &fld, sqNBASEm1);
    for (ib = 0; ib < sqNBASEm1; ib++)
      y2[ib] = base_from_char(fld[ib][0]);
    for (ib = 0; ib < sqNBASEm1; ib++) {
      read_next_line(f, buf);
      split_expecting(buf, &fld, sqNBASEm1+1);
      const base_t m = base_from_char(fld[0][0]);
      const base_t o = base_from_char(fld[0][1]);
      int jb;
      for (jb = 0; jb < sqNBASEm1; jb++) {
	const base_t n = y1[jb];
	const base_t p = y2[jb];
	const real_t val = value_from_string(fld[jb+1]);
	(*t)[i][j][k][l][m][n][o][p] = val;
      }
    }
  }
  fclose(f);
}

static void read_int21(int use_dna_params, int use_enthalpy_params, tab7_t *t)
{
  init_tab7(t,INF);
  FILE *f = parfile("int21", use_dna_params,use_enthalpy_params);
  char buf[MAXLINE+1];
  field_t fld;
  while (read_next_line(f, buf)) {
    if (!strstr(buf,"X"))
      die("read_int21: expecting 'X'");
    read_next_line(f, buf);
    if (split(buf, &fld) != 2)
      die("read_int21: expecting two bases");
    const base_t i = base_from_char(fld[0][0]);
    const base_t n = base_from_char(fld[1][0]);
    read_next_line(f, buf);
    if (split(buf, &fld) != 2)
      die("read_int21: expecting two bases");
    const base_t j = base_from_char(fld[0][0]);
    const base_t o = base_from_char(fld[1][0]);
    read_next_line(f, buf);
    if (split(buf, &fld) != 1)
      die("read_int21: expecting one field");
    if (fld[0][0] != 'Y')
      die("read_int21: expecting YB");
    const base_t m = base_from_char(fld[0][1]);
    read_next_line(f, buf);
    if (split(buf, &fld) != NBASE-1)
      die("read_int11: expecting %d bases", NBASE-1);
    base_t y[NBASE];
    int ib;
    for (ib = 0; ib < NBASE-1; ib++)
      y[ib] = base_from_char(fld[ib][0]);
    for (ib = 0; ib < NBASE-1; ib++) {
      read_next_line(f, buf);
      if (split(buf, &fld) != NBASE)
	die("read_int11: expecting %d fields", NBASE);
      const base_t k = base_from_char(fld[0][0]);
      int jb;
      for (jb = 0; jb < NBASE-1; jb++) {
	const base_t l = y[jb];
	const real_t val = value_from_string(fld[jb+1]);
	(*t)[i][j][k][l][m][n][o] = val;
      }
    }
  }
}

void param_read_from_text(const char *path, param_t p, 
                          int use_dna_params/*=0*/, int use_enthalpy_params/*=0*/)
{
  DIR *cwd = opendir(".");
  if (!cwd)
    die("param_read_from_text: cannot open current directory");
  if (chdir(path))
    die("param_read_from_text: cannot change directories to %s", path);
  p->use_dna_params = use_dna_params;
  p->use_enthalpy_params = use_enthalpy_params;
#define READ_STACK(x) read_stack(#x, p->use_dna_params, &p->x, use_enthalpy_params)
  READ_STACK(coaxial);
  READ_STACK(coaxstack);
  READ_STACK(stack);
  READ_STACK(tstack);
  READ_STACK(tstackcoax);
  READ_STACK(tstackh);
  READ_STACK(tstacki);
  READ_STACK(tstacki1n);
  READ_STACK(tstacki23);
  READ_STACK(tstackm);
#undef READ_STACK
  read_loop(p);
  read_miscloop(p);
  read_triloop(p);
  read_tloop(p);
  read_hexaloop(p);
  read_dangle(p->use_dna_params, p->use_enthalpy_params, &p->dangle_3p, &p->dangle_5p);
  read_int11(p->use_dna_params, p->use_enthalpy_params, &p->int11);
  read_int22(p->use_dna_params, p->use_enthalpy_params, &p->int22);
  read_int21(p->use_dna_params, p->use_enthalpy_params, &p->int21);
  if (fchdir(dirfd(cwd)))
    die("param_read_from_text: cannot cd back to original working directory");
  closedir(cwd);
}

//calculate gibbs free energy parameters at arbitrary temperature
//using enthalpy and gibbs free energy at 37C
//dG(T) = dH-T(dH=dG(37))/37
real_t dgT(real_t dg, real_t dh, real_t temperature){
  //convert dg and dh back to 37C
  if(isnan(dg)) return dg; //"infinite energy" is represented as nan
  double new_RT = R*temperature;
  dg = dg*RT;
  dh = dh*RT;
  real_t free_energy = dh - (dh-dg)*temperature/T37;
  return free_energy/new_RT;
}

void convert3(tab3_t p, tab3_t dg, tab3_t dh, real_t temperature,
                            int i, int j, int k)
{
  p[i][j][k] = dgT(dg[i][j][k],dh[i][j][k],temperature);
}

void convert4(tab4_t p, tab4_t dg, tab4_t dh, real_t temperature,
                            int i, int j, int k, int l)
{
  p[i][j][k][l] = dgT(dg[i][j][k][l],dh[i][j][k][l],temperature);
}

void param_read_alternative_temperature(const char *path, param_t p, 
                                        real_t temperature, int use_dna_params)
{
  //printf("Read free energy parameters");
  struct param DG;
  struct param DH;
  param_t dg=&DG;//gibbs free energy parameters at 37C
  param_t dh=&DH;//enthalpy parameters
  //printf("Read free energy parameters");
  param_read_from_text(path, p, use_dna_params,0);
  //printf("Read free energy parameters");
  param_read_from_text(path, &DG, use_dna_params,0);
  //printf("Read enthalpy parameters");
  param_read_from_text(path, &DH, use_dna_params,1);
  p->use_dna_params = use_dna_params;
  p->use_enthalpy_params = 0;
  p->temperature = temperature;
  int i,j,k,l,m,n,o,pp;
  for(i=0;i<NBASE;i++){
    for(j=0;j<NBASE;j++){
      for(k=0;k<NBASE;k++){
        convert3(p->dangle_3p,dg->dangle_3p,dh->dangle_3p,p->temperature,i,j,k);
        convert3(p->dangle_5p,dg->dangle_5p,dh->dangle_5p,p->temperature,i,j,k);
        for(l=0;l<NBASE;l++){
          convert4(p->coaxial,dg->coaxial,dh->coaxial,p->temperature,i,j,k,l);
          convert4(p->coaxstack,dg->coaxstack,dh->coaxstack,p->temperature,i,j,k,l);
          convert4(p->stack,dg->stack,dh->stack,p->temperature,i,j,k,l);
          convert4(p->tstack,dg->tstack,dh->tstack,p->temperature,i,j,k,l);
          convert4(p->tstackcoax,dg->tstackcoax,dh->tstackcoax,p->temperature,i,j,k,l);
          convert4(p->tstackh,dg->tstackh,dh->tstackh,p->temperature,i,j,k,l);
          convert4(p->tstacki,dg->tstacki,dh->tstacki,p->temperature,i,j,k,l);
          convert4(p->tstacki1n,dg->tstacki1n,dh->tstacki1n,p->temperature,i,j,k,l);
          convert4(p->tstacki23,dg->tstacki23,dh->tstacki23,p->temperature,i,j,k,l);
          convert4(p->tstackm,dg->tstackm,dh->tstackm,p->temperature,i,j,k,l);
          for(m=0;m<NBASE;m++){
            for(n=0;n<NBASE;n++){
              p->int11[i][j][k][l][m][n] = dgT(dg->int11[i][j][k][l][m][n],dh->int11[i][j][k][l][m][n],p->temperature);
              for(o=0;o<NBASE;o++){
                p->int21[i][j][k][l][m][n][o] = dgT(dg->int21[i][j][k][l][m][n][o],dh->int21[i][j][k][l][m][n][o],p->temperature);
                for(pp=0;pp<NBASE;pp++){
                  p->int22[i][j][k][l][m][n][o][pp] = dgT(dg->int22[i][j][k][l][m][n][o][pp],dh->int22[i][j][k][l][m][n][o][pp],p->temperature);
                }
              }
            }
          }
        }
      }
    }
  }
  for(i=0;i<p->ntriloop;i++)
    p->triloop[i].val = dgT(dg->triloop[i].val,dh->triloop[i].val,p->temperature);
  for(i=0;i<p->ntloop;i++)
    p->tloop[i].val = dgT(dg->tloop[i].val,dh->tloop[i].val,p->temperature);
  for(i=0;i<p->nhexaloop;i++)
    p->hexaloop[i].val = dgT(dg->hexaloop[i].val,dh->hexaloop[i].val,p->temperature);
  for(i=0;i<LOOP_MAX+1;i++){
    p->internal_loop_initiation[i] = dgT(dg->internal_loop_initiation[i],dh->internal_loop_initiation[i],p->temperature);
    p->bulge_loop_initiation[i] = dgT(dg->bulge_loop_initiation[i],dh->bulge_loop_initiation[i],p->temperature);
    p->hairpin_loop_initiation[i] = dgT(dg->hairpin_loop_initiation[i],dh->hairpin_loop_initiation[i],p->temperature);
  }
  p->Extrapolation_for_large_loops = dgT(dg->Extrapolation_for_large_loops,dh->Extrapolation_for_large_loops,p->temperature);
  p->maximum_correction = dgT(dg->maximum_correction,dh->maximum_correction,p->temperature);
  p->fm_array_first_element = dgT(dg->fm_array_first_element,dh->fm_array_first_element,p->temperature);
  p->a = dgT(dg->a,dh->a,p->temperature);
  p->b = dgT(dg->b,dh->b,p->temperature);
  p->c = dgT(dg->c,dh->c,p->temperature);
  p->multibranched_loop_offset = p->a;
  p->multibranched_loop_per_nuc_penalty = p->b;
  p->multibranched_loop_helix_penalty = p->c;
  p->a_2c = p->a + 2*p->c;
  p->a_2b_2c = p->a + 2*p->b + 2*p->c;
  p->terminal_AU_penalty = dgT(dg->terminal_AU_penalty,dh->terminal_AU_penalty,p->temperature);
  p->bonus_for_GGG_hairpin = dgT(dg->bonus_for_GGG_hairpin,dh->bonus_for_GGG_hairpin,p->temperature);
  p->c_hairpin_slope = dgT(dg->c_hairpin_slope,dh->c_hairpin_slope,p->temperature);
  p->c_hairpin_intercept = dgT(dg->c_hairpin_intercept,dh->c_hairpin_intercept,p->temperature);
  p->c_hairpin_of_3 = dgT(dg->c_hairpin_of_3,dh->c_hairpin_of_3,p->temperature);
  p->Bonus_for_Single_C_bulges_adjacent_to_C = dgT(dg->Bonus_for_Single_C_bulges_adjacent_to_C,dh->Bonus_for_Single_C_bulges_adjacent_to_C,p->temperature);

  p->Extrapolation_for_large_loops = dgT(dg->Extrapolation_for_large_loops,dh->Extrapolation_for_large_loops,p->temperature);
}



static void show1(const char *name, real_t t[], int n)
{
  int i;
  for (i = 0; i < n; i++)
    printf("%s[%d] = %.1f\n", name, i, RT*t[i]);
}

static void show3(const char *name, tab3_t * const t)
{
  int i[3];
  for (i[0] = 0; i[0] < NBASE; i[0]++)
    for (i[1] = 0; i[1] < NBASE; i[1]++)
      for (i[2] = 0; i[2] < NBASE; i[2]++)
	printf("%s[%c][%c][%c] = %.1f\n", 
	       name, 
	       base_as_char((base_t) i[0]), 
	       base_as_char((base_t) i[1]), 
	       base_as_char((base_t) i[2]),
	       RT*(*t)[i[0]][i[1]][i[2]]);
}

static void show4(const char *name, tab4_t * const t)
{
  int i[4];
  for (i[0] = 0; i[0] < NBASE; i[0]++)
    for (i[1] = 0; i[1] < NBASE; i[1]++)
      for (i[2] = 0; i[2] < NBASE; i[2]++)
	for (i[3] = 0; i[3] < NBASE; i[3]++)
	  printf("%s[%c][%c][%c][%c] = %.1f\n", 
		 name, 
		 base_as_char((base_t) i[0]), 
		 base_as_char((base_t) i[1]), 
		 base_as_char((base_t) i[2]),
		 base_as_char((base_t) i[3]),
		 RT*(*t)[i[0]][i[1]][i[2]][i[3]]);
}

static void show6(const char *name, tab6_t * const t)
{
  int i[6];
  for (i[0] = 0; i[0] < NBASE; i[0]++)
    for (i[1] = 0; i[1] < NBASE; i[1]++)
      for (i[2] = 0; i[2] < NBASE; i[2]++)
	for (i[3] = 0; i[3] < NBASE; i[3]++)
	  for (i[4] = 0; i[4] < NBASE; i[4]++)
	    for (i[5] = 0; i[5] < NBASE; i[5]++)
	      printf("%s[%c][%c][%c][%c][%c][%c] = %.1f\n", 
		     name, 
		     base_as_char((base_t) i[0]), 
		     base_as_char((base_t) i[1]), 
		     base_as_char((base_t) i[2]),
		     base_as_char((base_t) i[3]),
		     base_as_char((base_t) i[4]), 
		     base_as_char((base_t) i[5]), 
		     RT*(*t)[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]]);
}

static void show7(const char *name, tab7_t * const t)
{
  int i[7];
  for (i[0] = 0; i[0] < NBASE; i[0]++)
    for (i[1] = 0; i[1] < NBASE; i[1]++)
      for (i[2] = 0; i[2] < NBASE; i[2]++)
	for (i[3] = 0; i[3] < NBASE; i[3]++)
	  for (i[4] = 0; i[4] < NBASE; i[4]++)
	    for (i[5] = 0; i[5] < NBASE; i[5]++)
	      for (i[6] = 0; i[6] < NBASE; i[6]++)
		  printf("%s[%c][%c][%c][%c][%c][%c][%c] = %.1f\n", 
			 name, 
			 base_as_char((base_t) i[0]), 
			 base_as_char((base_t) i[1]), 
			 base_as_char((base_t) i[2]),
			 base_as_char((base_t) i[3]),
			 base_as_char((base_t) i[4]), 
			 base_as_char((base_t) i[5]), 
			 base_as_char((base_t) i[6]),
			 RT*(*t)[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]][i[6]]);
}

static void show8(const char *name, tab8_t * const t)
{
  int i[8];
  for (i[0] = 0; i[0] < NBASE; i[0]++)
    for (i[1] = 0; i[1] < NBASE; i[1]++)
      for (i[2] = 0; i[2] < NBASE; i[2]++)
	for (i[3] = 0; i[3] < NBASE; i[3]++)
	  for (i[4] = 0; i[4] < NBASE; i[4]++)
	    for (i[5] = 0; i[5] < NBASE; i[5]++)
	      for (i[6] = 0; i[6] < NBASE; i[6]++)
		for (i[7] = 0; i[7] < NBASE; i[7]++)
		  printf("%s[%c][%c][%c][%c][%c][%c][%c][%c] = %.1f\n", 
			 name, 
			 base_as_char((base_t) i[0]), 
			 base_as_char((base_t) i[1]), 
			 base_as_char((base_t) i[2]),
			 base_as_char((base_t) i[3]),
			 base_as_char((base_t) i[4]), 
			 base_as_char((base_t) i[5]), 
			 base_as_char((base_t) i[6]),
			 base_as_char((base_t) i[7]),
			 RT*(*t)[i[0]][i[1]][i[2]][i[3]][i[4]][i[5]][i[6]][i[7]]);
}

#define SHOW_SMALL_LOOP(X,N)			\
  static void show_##X##loop(const param_t p)	\
  {						\
    int i;					\
    printf("%s = %d\n", "n"#X"loop", p->n##X##loop);	\
    for (i = 0; i < p->n##X##loop; i++) {		\
      int j;						\
      for (j = 0; j < N; j++)				\
	printf("%s[%d].seq[%d] = %c\n", #X"loop", i, j, base_as_char(p->X##loop[i].seq[j])); \
      printf("%s[%d].val = %.1f\n", #X"loop", i, RT*p->X##loop[i].val);	\
    }									\
  }

SHOW_SMALL_LOOP(tri,5)
SHOW_SMALL_LOOP(t,6)
SHOW_SMALL_LOOP(hexa,8)

#undef SHOW_SMALL_LOOP

void param_show(const param_t p)
{
  printf("Using DNA parameters: %s\n", p->use_dna_params ? "yes" : "no");
#define SHOW(n,x) show##n(#x, &p->x)
  SHOW(4,coaxial);
  SHOW(4,coaxstack);
  SHOW(4,stack);
  SHOW(4,tstack);
  SHOW(4,tstackcoax);
  SHOW(4,tstackh);
  SHOW(4,tstacki);
  SHOW(4,tstacki1n);
  SHOW(4,tstacki23);
  SHOW(4,tstackm);
#define SHOW_LOOP(x) show1(#x, p->x, LOOP_MAX+1)
  SHOW_LOOP(internal_loop_initiation);
  SHOW_LOOP(bulge_loop_initiation);
  SHOW_LOOP(hairpin_loop_initiation);
#undef SHOW_LOOP
  printf("Extrapolation_for_large_loops = %g\n", RT*p->Extrapolation_for_large_loops);
#define SHOW_VAL(x) printf("%s = %.1f\n", #x, RT*p->x)
  SHOW_VAL(maximum_correction);
  SHOW_VAL(fm_array_first_element);
  SHOW_VAL(multibranched_loop_offset);
  SHOW_VAL(multibranched_loop_per_nuc_penalty);
  SHOW_VAL(multibranched_loop_helix_penalty);
  SHOW_VAL(terminal_AU_penalty);
  SHOW_VAL(bonus_for_GGG_hairpin);
  SHOW_VAL(c_hairpin_slope);
  SHOW_VAL(c_hairpin_intercept);
  SHOW_VAL(c_hairpin_of_3);
  SHOW_VAL(Bonus_for_Single_C_bulges_adjacent_to_C);
#undef SHOW_VAL
  show_triloop(p);
  show_tloop(p);
  show_hexaloop(p);
  SHOW(3,dangle_3p);
  SHOW(3,dangle_5p);
  SHOW(6,int11);
  SHOW(8,int22);
  SHOW(7,int21);
}

#define CHK 1.23456

void param_save_to_binary(const char *path, const param_t p)
{
  FILE *f = safe_fopen(path, "w");
  const real_t chk = CHK;
  fwrite(&chk, sizeof(chk), 1, f);
  fwrite(p, sizeof(struct param), 1, f);
  fclose(f);
}

void param_read_from_binary(const char *path, param_t p)
{
  FILE *f = safe_fopen(path, "r");
  real_t chk;
  fread(&chk, sizeof(chk), 1, f);
  if (chk != CHK)
    die("param_read_from_binary: %s: wrong format", path);
  fread(p, sizeof(struct param), 1, f);
  fclose(f);
}

#undef CHK


#ifdef OLD_FILE_FORMAT

static void read_stack_OLD(const char *name, int use_dna_params, tab4_t *p, int use_enthalpy_params)
{
  init_tab4(p);
  const int coaxial = !strcmp(name,"coaxial");
  FILE *f = parfile(name, use_dna_params, use_enthalpy_params);
  look_for_line_containing(f, "STACKING ");//used to stay "STACKING ENERGIES", altered to work with enthalpy
  char buf[MAXLINE+1];
  int i = 0, k = 0;
  while (fgets(buf, MAXLINE, f)) {
    if (strlen(buf) >= MAXLINE)
      die("read_stack: line too long");
    if (!strchr(buf,'.'))
      continue;
    int j = 0, l = 0;
    char *b;
    for (b = strtok(buf,whitespace); b; b = strtok(0,whitespace)) {
      if (coaxial)
	(*p)[j][i][k][l] = value_from_string(b);
      else
	(*p)[i][j][k][l] = value_from_string(b);
      l++;
      if (l == NBASE) {
	j++;
	l = 0;
      }
    }
    k++;
    if (k == NBASE) {
      i++;
      k = 0;
    }
  }
  fclose(f);
}

static void read_four_bases(const char *buf, base_t b[4])
{
  char a[4];
  const char *fmt = strchr(buf,'X') ? " %cX %cX %cX %cX" : " %c %c %c %c";
  if (sscanf(buf, fmt, &a[0], &a[1], &a[2], &a[3]) != 4)
    die("read_four_bases: error");
  int i;
  for (i = 0; i < 4; i++)
    b[i] = base_from_char(a[i]);
}

static void read_dangle_OLD(int use_dna_params, int use_enthalpy_params, tab3_t *d3p, tab3_t *d5p)
{
  init_tab3(d3p);
  init_tab3(d5p);
  FILE *f = parfile("dangle", use_dna_params, use_enthalpy_params);
  char buf[MAXLINE+1];
  while (fgets(buf, MAXLINE, f)) {
    if (strlen(buf) >= MAXLINE)
      die("read_dangle: line too long");
    if (is_only_whitespace(buf))
      continue;
    look_for_line_containing(f, "5' --> 3'");
    read_next_line(f, buf);
    tab3_t *d = strchr(buf,'X') ? d3p : d5p;
    base_t a[4], b[4];
    read_four_bases(buf, a);
    read_next_line(f, buf);
    read_four_bases(buf, b);
    expect_line_containing(f, "3' <-- 5'");
    real_t val[16];
    read_next_values(f, val, 16);
    int i, j;
    for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++)
	(*d)[a[i]][b[i]][j] = val[4*i+j];
  }
  fclose(f);
}

static void read_int11_OLD(int use_dna_params, int use_enthalpy_params, tab6_t *t)
{
  init_tab6(t);
  FILE *f = parfile("int11", use_dna_params, use_enthalpy_params);
  look_for_line_containing(f, "5' --> 3'");
  char buf[MAXLINE+1];
  while (fgets(buf, MAXLINE, f)) {
    if (strlen(buf) >= MAXLINE)
      die("read_int11: line too long");
    if (is_only_whitespace(buf))
      continue;
    look_for_line_containing(f, "5' --> 3'");
    expect_line_containing(f, "X");
    base_t b1[12], b2[12];
    read_next_line(f, buf);
    read_twelve_bases(buf, b1);
    read_next_line(f, buf);
    read_twelve_bases(buf, b2);
    expect_line_containing(f, "Y");
    expect_line_containing(f, "3' <-- 5'");
    int i;
    for (i = 0; i < 4; i++) {
      real_t val[24];
      read_next_values(f, val, 24);
      int j, k;
      for (j = 0; j < 6; j++)
	for (k = 0; k < 4; k++)
	  (*t)[b1[2*j]][i][b1[2*j+1]][b2[2*j+1]][k][b2[2*j]] = val[4*j+k];
    }
  }
  fclose(f);
}

static void read_six_bases(const char *buf, base_t b[6]) 
{
  char a[6];
  if (sscanf(buf, " Y%c Y%c Y%c Y%c Y%c Y%c", 
	     &a[0], &a[1], &a[2], &a[3], &a[4], &a[5]) != 6)
    die("read_six_bases: error");
  int i;
  for (i = 0; i < 6; i++)
    b[i] = base_from_char(a[i]);
}

static void read_twelve_bases(const char *buf, base_t b[12]) 
{
  char a[12];
  if (sscanf(buf, " %c %c %c %c %c %c %c %c %c %c %c %c", 
	     &a[0], &a[1], &a[2], &a[3],
	     &a[4], &a[5], &a[6], &a[7],
	     &a[8], &a[9], &a[10], &a[11]) != 12)
    die("read_twelve_bases: error");
  int i;
  for (i = 0; i < 12; i++)
    b[i] = base_from_char(a[i]);
}

static void read_int21_OLD(int use_dna_params, int use_enthalpy_params, tab7_t *t)
{
  init_tab7(t);
  FILE *f = parfile("int21", use_dna_params,use_enthalpy_params);
  look_for_line_containing(f, "5' --> 3'");
  char buf[MAXLINE+1];
  while (fgets(buf, MAXLINE, f)) {
    if (strlen(buf) >= MAXLINE)
      die("read_int21: line too long");
    if (is_only_whitespace(buf))
      continue;
    look_for_line_containing(f, "5' --> 3'");
    expect_line_containing(f, "X");
    base_t b1[12], b2[12], b3[6];
    read_next_line(f, buf);
    read_twelve_bases(buf, b1);
    read_next_line(f, buf);
    read_twelve_bases(buf, b2);
    read_next_line(f, buf);
    read_six_bases(buf, b3);
    expect_line_containing(f, "3' <-- 5'");
    int i;
    for (i = 0; i < 4; i++) {
      real_t val[24];
      read_next_values(f, val, 24);
      int j, k;
      for (j = 0; j < 6; j++)
	for (k = 0; k < 4; k++)
	  (*t)[b1[2*j]][b2[2*j]][i][k][b3[j]][b1[2*j+1]][b2[2*j+1]] = val[4*j+k];
    }
  }
  fclose(f);
}

static void read_int22_OLD(int use_dna_params, int use_enthalpy_params, tab8_t *t)
{
  init_tab8(t);
  FILE *f = parfile("int22", use_dna_params, use_enthalpy_params);
  look_for_line_containing(f, "5' ------> 3'");
  char buf[MAXLINE+1];
  while (fgets(buf, MAXLINE, f)) {
    if (strlen(buf) >= MAXLINE)
      die("read_int22: line too long");
    if (is_only_whitespace(buf))
      continue;
    look_for_line_containing(f, "5' ------> 3'");
    char a[4];
    if (!(fgets(buf, MAXLINE, f) && sscanf(buf, " %c \\/ \\_/ %c", &a[0], &a[1]) == 2))
      die("read_int22: couldn't read first line");
    if (!(fgets(buf, MAXLINE, f) && sscanf(buf, " %c /\\  |  %c", &a[2], &a[3]) == 2))
      die("read_int22: couldn't read second line");
    expect_line_containing(f, "3' <------ 5'");
    int i, j;
    base_t b[4];
    for (i = 0; i < 4; i++)
      b[i] = base_from_char(a[i]);
    for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++) {
	real_t val[16];
	read_next_values(f, val, 16);
	int k, l;
	for (k = 0; k < 4; k++)
	  for (l = 0; l < 4; l++)
	    (*t)[b[0]][b[1]][b[2]][b[3]][i][k][j][l] = val[4*k+l];
      }
  }
  fclose(f);
}

#endif /* OLD_FILE_FORMAT */


