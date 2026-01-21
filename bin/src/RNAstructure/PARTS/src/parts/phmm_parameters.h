#ifndef _PHMM_PARAMETERS_
#define _PHMM_PARAMETERS_

enum {STATE_ALN, STATE_INS1, STATE_INS2}; // States to access easier.

#define N_STATES (3)

static char state_names[][200] = {"STATE_ALN", "STATE_INS1", "STATE_INS2"};

#endif // _PHMM_PARAMETERS_
