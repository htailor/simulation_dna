#ifndef INCLUDE_MATRIX_HH
#define INCLUDE_MATRIX_HH

#include <gsl/gsl_matrix.h>

class Matrix{
 
public:

 Matrix();
 Matrix(int n_row, int n_col);
 int n_row();
 int n_column();
 double get(int i, int j);
 void set(int i, int j, double value);
 void print_to_file(const char* filename);
 void print_row(const char* filename, int ROW);
 void print_to_file_list(const char* filename);
 gsl_matrix return_gsl_matrix();


private:

 int _N_ROW;
 int _N_COLUMN;

 gsl_matrix *_MATRIX;

};

#endif
