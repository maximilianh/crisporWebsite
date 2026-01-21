#include <stdio.h>
#include <stdlib.h>
#include "phmm_array.h"
#include "phmm.h"
#include "structure/structure_object.h"
#include "utils/xmath/log/xlog_math.h"

bool _DUMP_PHMM_ARRAY_MESSAGES_ = false;

t_phmm_array::t_phmm_array(int _n1, int _n2, int _phmm_band_constraint_size, bool mallocate)
{
	n1 = _n1;
	n2 = _n2;

if(_DUMP_PHMM_ARRAY_MESSAGES_)
	printf("Allocing phmm array..\n");

	this->n_bytes_alloced = 0.0f;

	if(mallocate)
	{
		this->array = (double***)malloc(sizeof(double**) * (n1 + 2));
	}
	else
	{
		this->array = NULL;
	}

	this->phmm_band_constraint_size = _phmm_band_constraint_size;
	this->set_hmm_array_banded_limits();

	n_bytes_alloced += (sizeof(double**) * (n1 + 2));

	//for(int i = 0; i <= n1 + 1; i++)
	//{
	//	//printf("%d -> %d, %d\n", i, this->low_phmm_array_limits[i], this->high_phmm_array_limits[i]);
	//}
	//

	for(int i = 0; i <= n1 + 1; i++)
	{
		int low_k = this->low_phmm_array_limits[i];
		int high_k = this->high_phmm_array_limits[i];
		//printf("%d -> %d, %d\n", i, low_k, high_k);
		//int low_k = 0;
		//int high_k = n2+1;

		if(mallocate)
		{
			this->array[i] = (double**)malloc(sizeof(double*) * (n2 + 2));
			this->array[i] -= low_k;
		}

		this->n_bytes_alloced += (sizeof(double*) * (high_k - low_k + 1));		

if(_DUMP_PHMM_ARRAY_MESSAGES_)
		printf("At %lf bytes for phmm array.\r", n_bytes_alloced);

		for(int k = low_k; k <= high_k; k++)
		{
			if(mallocate)
			{
				this->array[i][k] = (double*)malloc(sizeof(double) * (N_STATES + 5));
				//printf("%d, %d -> %d\n", i, k, this->array[i][k]);
			}

			this->n_bytes_alloced += (sizeof(double) * N_STATES);

			for(int i_state = STATE_INS1; mallocate && i_state <= STATE_ALN; i_state++)
			{
				this->array[i][k][i_state] = xlog(0.0f);
			}
		} // k loop
	} // i loop

	//this->array[2][6][0] = 1221233214.0f;
	//if(this->array[2][8][0] == 1221233214.0f)
	//{
	//	printf("WTF1! (%d, %d) \n", (int)this->array[2][6] , (int)this->array[2][8]);
	//	exit(0);
	//}

	//getc(stdin);

if(_DUMP_PHMM_ARRAY_MESSAGES_)
	printf("%lf bytes allocated for phmm_array\n", this->n_bytes_alloced);
}

t_phmm_array::~t_phmm_array()
{
	if(this->array != NULL)
	{
		//this->array = (double***)malloc(sizeof(double**) * (n1 + 4));
		for(int i = 0; i <= n1 + 1; i++)
		{
			//this->array[k] = (double**)malloc(sizeof(double*) * (n2 + 4));
			int low_k = this->low_phmm_array_limits[i];
			int high_k = this->high_phmm_array_limits[i];

			for(int k = low_k; k <= high_k; k++)
			{
				//this->array[k][k] = (double**)malloc(sizeof(double) * N_STATES);
				free(this->array[i][k]);
			}

			this->array[i] += low_k;
			free(this->array[i]);
		}

		free(this->array);
	}

	free(this->low_phmm_array_limits);
	free(this->high_phmm_array_limits);
}

bool t_phmm_array::check_phmm_boundary(int i, int k)
{
	if(this->low_phmm_array_limits[i] <= k &&
		this->high_phmm_array_limits[i] >= k)
	{
		return(true);
	}

	return(false);
}

int t_phmm_array::low_phmm_limit(int i, int n1, int n2, int phmm_band_constraint_size)
{
	if(i == n1+1)
	{
		return(n2+1);
	}

	int corresponding_i_in_seq2 = (int)(((double)i * (double)n2) / (double)n1);
	return(MAX(0, corresponding_i_in_seq2 - phmm_band_constraint_size));
}

int t_phmm_array::high_phmm_limit(int i, int n1, int n2, int phmm_band_constraint_size)
{
	if(i == n1+1)
	{
		return(n2+1);
	}

	int corresponding_i_in_seq2 = (int)(((double)i * (double)n2) / (double)n1);
	return(MIN(n2, corresponding_i_in_seq2 + phmm_band_constraint_size));
}

bool t_phmm_array::check_phmm_boundary(int i, int k, int n1, int n2, int phmm_band_constraint_size)
{
	if(k <= t_phmm_array::high_phmm_limit(i, n1, n2, phmm_band_constraint_size) &&
		k >= t_phmm_array::low_phmm_limit(i, n1, n2, phmm_band_constraint_size))
	{
		return(true);
	}
	
	return(false);
}

void t_phmm_array::set_hmm_array_banded_limits()
{
	this->low_phmm_array_limits = (int*)malloc(sizeof(int) * (this->n1 + 2));
	this->high_phmm_array_limits = (int*)malloc(sizeof(int) * (this->n1 + 2));

	for(int i = 0; i <= this->n1+1; i++)
	{
		this->low_phmm_array_limits[i] = t_phmm_array::low_phmm_limit(i, n1, n2, this->phmm_band_constraint_size);
		this->high_phmm_array_limits[i] = t_phmm_array::high_phmm_limit(i, n1, n2, this->phmm_band_constraint_size);

if(_DUMP_PHMM_ARRAY_MESSAGES_)
		printf("%d: %d-%d\r", i, this->low_phmm_array_limits[i], this->high_phmm_array_limits[i]);
	} // i loop
}


double& t_phmm_array::x(int i, int k, int state)
{
	return(this->array[i][k][state]);
}


