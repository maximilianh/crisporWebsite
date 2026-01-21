#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "param.h"
#include "util.h"

static const char Usage[] = 
  "Usage: %s <file>\n"
  " or    %s -v\n"
  "\n"
  "Writes parameters from the directory given by $DATAPATH to <file>,\n"
  "in native binary format\n"
  "\n"
  "options:\n"
  "-d: use DNA parameters\n"
  "-v: just show parameters\n\n";

int main(int argc, char **argv)
{
  const char *cmd = *argv;
  int use_dna_params = 0, verbose = 0;

  /* process command-line arguments */
  int c;
  while ((c = getopt(argc, argv, "hdv")) != EOF)
    if (c == 'h')
      die(Usage,cmd,cmd);
    else if (c == 'd')
      use_dna_params = 1;
    else if (c == 'v')
      verbose = 1;
    else
      die(Usage,cmd,cmd);
  argc -= optind;
  argv += optind;

  if (!verbose && argc < 1)
    die(Usage,cmd,cmd);

  /* read parameters */
  const char *path = getenv("DATAPATH");
  if (!path)
    die("%s: need to set environment variable $DATAPATH", cmd);
  struct param p;
  param_read_from_text(path, &p, use_dna_params, 0);
  
  /* write parameters */
  if (verbose)
    param_show(&p);
  else
    param_save_to_binary(*argv, &p);
  
  return 0;
  
}
