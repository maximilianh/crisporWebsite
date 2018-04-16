extern "C"
{

/* functions from part_func.c */
float  pf_fold(char *sequence, char *structure);
/* calculate partition function and base pair probabilities */
void   init_pf_fold(int length);    /* allocate space for pf_fold() */
void   free_pf_arrays(void);        /* free arrays from pf_fold() */
void   update_pf_params(int length); /*recalculate energy parameters */
char   bppm_symbol(float *x);   /* string representation of structure */
double mean_bp_dist(int length); /* mean pair distance of ensemble */

}
