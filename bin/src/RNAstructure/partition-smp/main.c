#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <sys/stat.h>
#include <dirent.h>
#include <ctype.h>

#include "prna.h"
#include "util.h"
#include "base.h"


static const char Usage[] = 
  "Usage: %s [options] <sequence, or file containing one>\n"
  "\n"
  "options:\n"
  "-h:        show this message\n"
  "-b <file>: read parameters from <file>, in native binary format\n"
  "-d:        use DNA parameters\n"
  "-t <file>: write probability matrix as text to <file>\n"
  "-l <file>: write -log10 probabilities as text to <file>\n"
  "-p <file>: write ProbKnot structure in ct format to <file>\n"
  "-m <length>: set minimum helix length for ProbKnot\n"
  "             (default: 3 base pairs)\n"
  "-v:        show arrays\n\n"
  "If none of -t, -l, -p, or -v is chosen,\n"
  "writes ProbKnot structure in ct format to stdout\n";

inline static int cp(int i, int j, const base_t *s)
{
  return j-i-1 >= LOOP_MIN && is_canonical_pair(s[i],s[j]);
}

inline static int can_pair(int i, int j, int n, const base_t *s, const char *seq)
{
  if (!isupper(seq[i]) || !isupper(seq[j]))
    return 0;
  if (j < i) {
    const int tmp = i;
    i = j;
    j = tmp;
  }
  return cp(i,j,s) && ((i > 0 && j < n-1 && cp(i-1,j+1,s)) || cp(i+1,j-1,s));
}  
/*
void test_bcp(char *s){
  int n = strlen(s);
  
  int i, j;
  //int *neighboring_pair;

  base_t *seq = (base_t *) safe_malloc(n*sizeof(base_t));
  sequence_from_string(seq, s);

  for (i=0; i<n; i++){
    for (j=i+1; j<n; j++){
      if ((j-i < LOOP_MIN+1) || !isupper(s[i]) || !isupper(s[j])){
        if (can_pair(i,j,n,seq,s)!=0)
          printf("%d (%c) : %d (%c) Shouldn't Pair\n", i, s[i], j, s[j]);
      }
      else{

        if (can_pair(i,j,n,seq,s) != is_canonical_pair(seq[i],seq[j]) && ((i > 0 && j < n-1 && is_canonical_pair(seq[i-1],seq[j+1])) || (j-i>=LOOP_MIN+3 && is_canonical_pair(seq[i+1],seq[j-1]))))
          printf("%d (%c) : %d (%c) can_pair( %d ), bcp( %d )\n", i, s[i], j, s[j], can_pair(i,j,n,seq,s), is_canonical_pair(seq[i],seq[j]) && ((i > 0 && j < n-1 && is_canonical_pair(seq[i-1],seq[j+1])) || (is_canonical_pair(seq[i+1],seq[j-1]))));
      }
    }
  }
} 
//*/
int main(int argc, char **argv)
{
  const char *cmd = *argv;
  int use_dna_params = 0, min_helix_length = 3, verbose = 0;
  const char *neg_log10_filename = 0;
  const char *text_matrix_filename = 0;
  const char *probknot_filename = 0; 
  const char *binary_parameter_filename = 0;

  /* process command-line arguments */
  int c;
  while ((c = getopt(argc, argv, "hb:dt:l:p:m:v")) != EOF)
    if (c == 'h')
      die(Usage,cmd);
    else if (c == 'b')
      binary_parameter_filename = optarg;
    else if (c == 'd')
      use_dna_params = 1;
    else if (c == 't')
      text_matrix_filename = optarg;
    else if (c == 'l')
      neg_log10_filename = optarg;
    else if (c == 'p')
      probknot_filename = optarg;
    else if (c == 'm')
      min_helix_length = atoi(optarg);
    else if (c == 'v')
      verbose = 1;
    else
      die(Usage,cmd);
  argc -= optind;
  argv += optind;
  if (argc == 0)
    die(Usage,cmd);

  /* get sequence */
  char *seq = sequence(*argv);
  
  /* read parameters */
  struct param par;
  if (binary_parameter_filename) {
    param_read_from_binary(binary_parameter_filename, &par);
    if (par.use_dna_params != use_dna_params)
      die("%s: -d option %s, but '%s' %s DNA parameters", cmd,
	  use_dna_params ? "set" : "not set", 
	  binary_parameter_filename, 
	  par.use_dna_params ? "is from" : "is not from");
  } else {
    const char *path = getenv("DATAPATH");
    if (!path)
      die("%s: need to set environment variable $DATAPATH", cmd);
    param_read_from_text(path, &par, use_dna_params,0);
  }
  
  /* calculate partition function */
  //test_bcp(seq);
  
  prna_t p = prna_new(seq, &par,!verbose, generate_bcp(seq));
  
  //p->base_can_pair=generate_bcp(seq);

  /* output */
  if (neg_log10_filename){
  //  printf("Writing Probabilities");
    prna_write_neg_log10_probabilities(p,neg_log10_filename);
  //  printf("DoneWriting Probabilities");
  }
  if (text_matrix_filename)
    prna_write_probability_matrix(p,text_matrix_filename);
  if (probknot_filename)
    prna_write_probknot(p,probknot_filename,seq,min_helix_length);
  if (verbose)
    prna_show(p);

  if (!(neg_log10_filename ||
	text_matrix_filename ||
	probknot_filename || 
	verbose))
    prna_write_probknot(p,0,seq,min_helix_length);

  /* cleanup */
  //printf("Cleaning Up");
  prna_delete(p);
  free(seq);
  
  return 0;
  
}
