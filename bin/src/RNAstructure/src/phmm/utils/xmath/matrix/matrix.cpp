#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "../../file/utils.h"

bool _DUMP_MATRIX_MESSAGES_ = false;

t_matrix::t_matrix(int _height, int _width, bool _symmetric)
{
	this->symmetric = _symmetric;
	this->height = _height;
	this->width = _width;

	// Allocate array.
	this->allocate_matrix(NULL);
}

// Default constructor.
t_matrix::t_matrix()
{
	this->symmetric = false;
	this->height = 0;
	this->width = 0;
	this->matrix = NULL;
}

t_matrix::t_matrix(double** _matrix, int _height, int _width, bool _symmetric)
{
	this->symmetric = _symmetric;
	this->height = _height;
	this->width = _width;

	// Allocate array.
	this->allocate_matrix(_matrix);
}

// Copy constructor.
t_matrix::t_matrix(t_matrix* _matrix)
{
	this->symmetric = _matrix->symmetric;
	this->height = _matrix->height;
	this->width = _matrix->width;

	// Allocate array.
	this->allocate_matrix(_matrix->matrix);
}

t_matrix::~t_matrix()
{
	// Free array.
	//this->matrix = (double**)malloc(sizeof(double*) * (this->height + 2));
	for(int i_row = 0; i_row <= this->height; i_row++)
	{		
		if(this->symmetric)
		{
			// Note that there is no pointer shift here.
			//this->matrix[i_row] = (double*)malloc(sizeof(double) * (this->width + 2));

			// Pointer shift.
			//this->matrix[i_row] -= i_row;
			this->matrix[i_row] += i_row;
			free(this->matrix[i_row]);
		}
		else
		{
			//this->matrix[i_row] = (double*)malloc(sizeof(double) * (this->width + 2));			
			free(this->matrix[i_row]);
		} // symmetry check.
	} // i_row loop

	free(this->matrix);
}

void t_matrix::init_by_constant(double init_value)
{
	for(int i_row = 1; i_row <= this->height; i_row++)
	{		
		if(this->symmetric)
		{
			for(int i_col = i_row; i_col <= this->width; i_col++)			
			{
				this->x(i_row, i_col) = init_value;
			} // i_col loop
		}
		else
		{
			for(int i_col = 1; i_col <= this->width; i_col++)			
			{
				this->x(i_row, i_col) = init_value;
			} // i_col loop
		}
	} // i_row loop
}

void t_matrix::add(t_matrix* _m2add)
{
	for(int i_row = 1; i_row <= this->height; i_row++)
	{		
		if(this->symmetric)
		{
			for(int i_col = i_row; i_col <= this->width; i_col++)			
			{
				this->x(i_row, i_col) += _m2add->x(i_row, i_col);
			} // i_col loop
		}
		else
		{
			for(int i_col = 1; i_col <= this->width; i_col++)			
			{
				this->x(i_row, i_col) += _m2add->x(i_row, i_col);
			} // i_col loop
		}
	} // i_row loop
}

void t_matrix::sub(t_matrix* _m2sub)
{
	for(int i_row = 1; i_row <= this->height; i_row++)
	{		
		if(this->symmetric)
		{
			for(int i_col = i_row; i_col <= this->width; i_col++)			
			{
				this->x(i_row, i_col) -= _m2sub->x(i_row, i_col);
			} // i_col loop
		}
		else
		{
			for(int i_col = 1; i_col <= this->width; i_col++)			
			{
				this->x(i_row, i_col) -= _m2sub->x(i_row, i_col);
			} // i_col loop
		}
	} // i_row loop
}

void t_matrix::mul(double val)
{
	for(int i_row = 1; i_row <= this->height; i_row++)
	{		
		if(this->symmetric)
		{
			for(int i_col = i_row; i_col <= this->width; i_col++)			
			{
				this->x(i_row, i_col) *= val;
			} // i_col loop
		}
		else
		{
			for(int i_col = 1; i_col <= this->width; i_col++)			
			{
				this->x(i_row, i_col) *= val;
			} // i_col loop
		}
	} // i_row loop
}

/*
Replace 0's by minimum value of the matrix.
*/
void t_matrix::fix_zeros_by_eps(double zero_epsilon)
{
	for(int i_row = 1; i_row <= this->height; i_row++)
	{		
		if(this->symmetric)
		{
			for(int i_col = i_row; i_col <= this->width; i_col++)			
			{
				if(this->x(i_row, i_col) < zero_epsilon)
				{
					this->x(i_row, i_col) = zero_epsilon;
				}
			} // i_col loop
		}
		else
		{
			for(int i_col = 1; i_col <= this->width; i_col++)			
			{
				if(this->x(i_row, i_col) < zero_epsilon)
				{
					this->x(i_row, i_col) = zero_epsilon;
				}
			} // i_col loop
		}
	} // i_row loop
}

double& t_matrix::x(int i_row, int i_col)
{
	//if(i_row >= 1 && i_row <= this->height &&
	//	i_col >= 1 && i_col <= this->width)
	{
		if(this->symmetric)
		{
			if(i_col >= i_row)
			{
				return(this->matrix[i_row][i_col]);
			}
			else
			{
				return(this->matrix[i_col][i_row]);
			}			
		} // symmetric?
		else
		{
			return(this->matrix[i_row][i_col]);
		} // not symmetric?
	}
	//else
	//{
	//	int* p = NULL;
	//	*p = 0;
	//	printf("Accessing indices out of matrix size.\n");
	//	exit(0);
	//}
}

void t_matrix::dump_matrix(char* fp)
{
	FILE* f_matrix = open_f(fp, "w");

	printf("Dumping to %s\n", fp);

	// Dump indices are 1-based.
	for(int i_row = 1; i_row <= this->height; i_row++)
	{
		for(int i_col = 1; i_col <= this->width; i_col++)
		{
			fprintf(f_matrix, "%lf ", this->x(i_row, i_col));
		} // i_col loop

		fprintf(f_matrix, "\n");
	} // i_row loop

	fclose(f_matrix);
}

void t_matrix::dump_sparse_matrix(char* fp)
{
	FILE* f_matrix = open_f(fp, "wb");

	// Must dump all the entries without regard to symmetry of the matrix.
	for(int i_row = 1; i_row <= this->height; i_row++)
	{
		for(int i_col = 1; i_col <= this->width; i_col++)
		{
			if(i_row > i_col && this->symmetric)
			{
				double cur_val = this->x(i_col, i_row);
				fwrite((void*)&i_row, sizeof(int), 1, f_matrix);
				fwrite((void*)&i_col, sizeof(int), 1, f_matrix);
				fwrite((void*)&cur_val, sizeof(double), 1, f_matrix);
			}
			else
			{
				double cur_val = this->x(i_row, i_col);
				fwrite((void*)&i_row, sizeof(int), 1, f_matrix);
				fwrite((void*)&i_col, sizeof(int), 1, f_matrix);
				fwrite((void*)&cur_val, sizeof(double), 1, f_matrix);
			}
			
		} // i_col loop
	} // i_row loop

	fclose(f_matrix);
}

void t_matrix::load_sparse_matrix(char* fp)
{
	FILE* f_matrix = open_f(fp, "rb");

	int cur_i;
	int cur_j; 
	double cur_value;
	//while(fscanf(f_matrix, "%d %d %lf", &cur_i, &cur_j, &cur_value) == 3)
	while(fread(&cur_i, sizeof(int), 1, f_matrix) == 1)
	{
		if(fread(&cur_j, sizeof(int), 1, f_matrix) != 1)
		{
			printf("Could not read current j in %s @ %s(%d)\n", fp, __FILE__, __LINE__);
			exit(0);
		}

		if(fread(&cur_value, sizeof(double), 1, f_matrix) != 1)
		{
			printf("Could not read current value in %s @ %s(%d)\n", fp, __FILE__, __LINE__);
			exit(0);
		}

		//printf("Read %d, %d %lf\n", cur_i, cur_j, cur_value);

		// If the matrix is symmetric, do a check on the read indices.
		if(this->symmetric)
		{
			if(cur_j > cur_i)
			{
				this->x(cur_i, cur_j) = cur_value;
			}
		}
		else
		{
			this->x(cur_i, cur_j) = cur_value;
		}
	} // file reading loop.

	fclose(f_matrix);
}

double t_matrix::correlate(double** a_matrix)
{
	double correlation = 0.0f;
	for(int i = 1; i <= this->height; i++)
	{
		for(int j = 1; j <= this->width; j++)
		{
			correlation += this->x(i,j) * a_matrix[i][j];
		} // j loop
	} // i loop

	return(correlation);
}

double t_matrix::correlate(t_matrix* a_matrix)
{
	double correlation = 0.0f;
	for(int i = 1; i <= this->height; i++)
	{
		for(int j = 1; j <= this->width; j++)
		{
			correlation += this->x(i,j) * a_matrix->x(i,j);
		} // j loop
	} // i loop

	return(correlation);
}

t_matrix* t_matrix::correlation_matrix(t_matrix* a_matrix)
{
	t_matrix* correlation_matrix = new t_matrix(this->height, this->width, false);
	for(int i = 1; i <= this->height; i++)
	{
		for(int j = 1; j <= this->width; j++)
		{
			correlation_matrix->x(i, j) = this->x(i, j) * a_matrix->x(i, j);
		} // j loop
	} // i loop

	return(correlation_matrix);
}

t_matrix* t_matrix::correlation_matrix(double** a_matrix)
{
	t_matrix* correlation_matrix = new t_matrix(this->height, this->width, false);
	for(int i = 1; i <= this->height; i++)
	{
		for(int j = 1; j <= this->width; j++)
		{
			correlation_matrix->x(i, j) = this->x(i, j) * a_matrix[i][j];
		} // j loop
	} // i loop

	return(correlation_matrix);
}

void t_matrix::allocate_matrix(double** initing_matrix)
{
	// Set memory usage to 0 bytes.
	this->mem_usage = 0.0f;

	if(this->symmetric && this->width != this->height)
	{
		printf("Cannot allocate a symmetric matric with unequal width and height\n");
		exit(0);
	}

	this->matrix = (double**)malloc(sizeof(double*) * (this->height + 2));
	this->mem_usage += (sizeof(double*) * (this->height + 2));

	// Allocate matrix: Make sure allocation starts from 0 since free'ing starts at 0.
	for(int i_row = 0; i_row <= this->height; i_row++)
	{		
		if(this->symmetric)
		{
			// Use the pointer shift and do not allocate symmetric indices.
			this->matrix[i_row] = (double*)malloc(sizeof(double) * (this->width - i_row + 2));
			this->mem_usage += (sizeof(double) * (this->width - i_row + 2)); // Save memory for symmetric matrices.
			//this->mem_usage += (sizeof(double) * (this->width + 2));

			// Pointer shift.
			this->matrix[i_row] -= i_row;

			for(int i_col = i_row; i_col <= this->width; i_col++)			
			{
				this->matrix[i_row][i_col] = 0.0f;
			} // i_col loop
		} // symmetry check.
		else
		{
			this->matrix[i_row] = (double*)malloc(sizeof(double) * (this->width + 2));
			this->mem_usage += (sizeof(double) * (this->width + 2));

			for(int i_col = 0; i_col <= this->width; i_col++)
			{
				this->matrix[i_row][i_col] = 0.0f;
			} // i_col loop
		}
	} // i_row loop

	// Init. matrix.
	for(int i_row = 0; i_row <= this->height; i_row++)
	{		
		if(this->symmetric)
		{
			for(int i_col = i_row; i_col <= this->width; i_col++)			
			{
				if(initing_matrix != NULL)
				{
					this->matrix[i_row][i_col] = initing_matrix[i_row][i_col];
				}
				else
				{
					this->matrix[i_row][i_col] = 0.0f;
				}
			} // i_col loop
		} // symmetry check.
		else
		{
			for(int i_col = 0; i_col <= this->width; i_col++)
			{
				if(initing_matrix != NULL)
				{
					this->matrix[i_row][i_col] = initing_matrix[i_row][i_col];
				}
				else
				{
					this->matrix[i_row][i_col] = 0.0f;
				}
			} // i_col loop
		}
	} // i_row loop
}

double** t_matrix::get_matrix()
{
	return(this->matrix);
}

void t_matrix::normalize_by_max()
{
	double matrix_max = -1000000.0f;
	for(int i_row = 1; i_row <= this->height; i_row++)
	{		
		for(int i_col = 1; i_col <= this->width; i_col++)			
		{
			//printf("(%d, %d) = %.10f\n", i_row, i_col, this->x(i_row, i_col));
			if(this->x(i_row, i_col) > matrix_max)
			{
				matrix_max = this->x(i_row, i_col);
			}
		} // i_col loop.
	} // i_row loop.

	// Matrix maximum is determined, normalize: If the matrix maximum is identically 0.0, the 
	// normalization cannot be done, in this case 0/0 is treated as 1 and all matrix entries are set to 1.0.
	if(matrix_max == 0.0f)
	{
		for(int i_row = 1; i_row <= this->height; i_row++)
		{	
			if(this->symmetric)
			{
				for(int i_col = i_row; i_col <= this->width; i_col++)			
				{
					this->x(i_row, i_col) = 1.0f;
				} // i_col loop.

			}
			else
			{
				for(int i_col = 1; i_col <= this->width; i_col++)			
				{
					this->x(i_row, i_col) = 1.0f;
				} // i_col loop.
			}
		} // i_row loop.
	} // check if matrix_max is identically 0.0.
	else
	{
if(_DUMP_MATRIX_MESSAGES_)
		printf("Matrix max is %.10f\n", matrix_max);
		//getc(stdin);
		for(int i_row = 1; i_row <= this->height; i_row++)
		{	
			if(this->symmetric)
			{
				for(int i_col = i_row; i_col <= this->width; i_col++)			
				{
					this->x(i_row, i_col) /= matrix_max;
				} // i_col loop.

			}
			else
			{
				for(int i_col = 1; i_col <= this->width; i_col++)			
				{
					this->x(i_row, i_col) /= matrix_max;
				} // i_col loop.
			}
		} // i_row loop.
	} // check if matrix_max is not identically 0.0.
}

void t_matrix::powerize_each_element(double power_factor)
{
	for(int i_row = 1; i_row <= this->height; i_row++)
	{		
		if(this->symmetric)
		{
			for(int i_col = i_row; i_col <= this->width; i_col++)			
			{
				this->x(i_row, i_col) = pow(this->x(i_row, i_col), power_factor);
			} // i_col loop.
		}
		else
		{
			for(int i_col = 1; i_col <= this->width; i_col++)			
			{
				this->x(i_row, i_col) = pow(this->x(i_row, i_col), power_factor);
			} // i_col loop.
		}
	} // i_row loop.
}

void t_matrix::exponentiate_by_element(double base)
{
	for(int i_row = 1; i_row <= this->height; i_row++)
	{		
		if(this->symmetric)
		{
			for(int i_col = i_row; i_col <= this->width; i_col++)			
			{
				//this->x(i_row, i_col) = exp(this->x(i_row, i_col));
				this->x(i_row, i_col) = pow(base, this->x(i_row, i_col));
			} // i_col loop.
		}
		else
		{
			for(int i_col = 1; i_col <= this->width; i_col++)			
			{
				//this->x(i_row, i_col) = exp(this->x(i_row, i_col));
				this->x(i_row, i_col) = pow(base, this->x(i_row, i_col));
			} // i_col loop.
		}
	} // i_row loop.
}