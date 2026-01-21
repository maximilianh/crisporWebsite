#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <sys/stat.h>
#include <dirent.h>

#include "frna.h"
#include "fbase.h"
#include "../partition-smp/util.h"

static const char Usage[] =
  "Usage: %s [options] <sequence, or file containing one> <output file>\n"
  "\n"
  "options:\n"
  "-h:        show this message\n"
  "-d:        use DNA parameters\n"
  "-v:        show arrays\n\n"
  "the output file is a fold save file\n"
  "to perform traceback for structures, run the refold program on this file";

int main(int argc, char **argv)
{
  const char *cmd = *argv;
  int use_dna_fparams = 0, min_helix_length = 3, verbose = 0;
  /* process command-line arguments */
  int c;
  while ((c = getopt(argc, argv, "h:d:v")) != EOF)
    if (c == 'h')
      die(Usage,cmd);
    else if (c == 'd')
      use_dna_fparams = 1;
    else if (c == 'v')
      verbose = 1;
    else
      die(Usage,cmd);
  argc -= optind;
  argv += optind;
  if (argc != 2)
    die(Usage,cmd);
  /* get sequence */
  char *seq = fsequence(*argv);
  char *outfile = argv[1];
  /* read fparameters */
  struct fparam par;
  
  const char *path = getenv("DATAPATH");
  if (!path)
    die("%s: need to set environment variable $DATAPATH", cmd);
  fparam_read_from_text(path, use_dna_fparams, &par);
 
// fparam_show(par);
  /* calculate partition function */
  frna_t p = frna_new(seq, &par);
  if(verbose) frna_show(p);
  frna_write_save_file(p,&par,outfile);

  /* cleanup */
  frna_delete(p);
  free(seq);

  return 0;

}
