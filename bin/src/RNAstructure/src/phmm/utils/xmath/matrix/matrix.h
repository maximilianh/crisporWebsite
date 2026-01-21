#ifndef __MATRIX__
#define __MATRIX__

#include <math.h>

/*
Accession similar to MATLAB matrices.
Start all indices from 1. The operations, normalizing and powerizing are done with 1 based indexing.
*/
class t_matrix
{
public:
	t_matrix(int _height, int _width, bool _symmetric);
	t_matrix(double** _matrix, int _height, int _width, bool _symmetric);
	t_matrix(t_matrix* _matrix);
	t_matrix();
	~t_matrix();

	int height;
	int width;
	bool symmetric;
	double mem_usage;

	void allocate_matrix(double** matrix);

	// Accessor/Setter
	double& x(int i_row, int i_col);

	double** get_matrix();

	void fix_zeros_by_eps(double zero_eps);

	void dump_matrix(char* fp);
	void dump_sparse_matrix(char* fp);
	void load_sparse_matrix(char* fp);

	void normalize_by_max();
	void powerize_each_element(double power_factor);
	void add(t_matrix* _m2add);
	void sub(t_matrix* _m2sub);
	void mul(double val);
	void init_by_constant(double init_value);

	// Take exp of each element.
	void exponentiate_by_element(double base);

	double correlate(double** a_matrix);
	double correlate(t_matrix* a_matrix);
	t_matrix* correlation_matrix(t_matrix* a_matrix);
	t_matrix* correlation_matrix(double** a_matrix);

private:
	double** matrix;
};

#endif // __MATRIX__
