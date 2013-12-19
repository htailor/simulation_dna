#include "Matrix.hh"


using namespace std;

Matrix::Matrix(){}

Matrix::Matrix(int n_row, int n_column){
   
   _N_ROW = n_row;
   _N_COLUMN = n_column;

   _MATRIX = gsl_matrix_alloc(_N_ROW,_N_COLUMN);
}

int Matrix::n_row(){
   return _N_ROW;
}

int Matrix::n_column(){
   return _N_COLUMN;
}

void Matrix::set(int i,int j, double value){
      gsl_matrix_set(_MATRIX,i,j,value); 
}

double Matrix::get(int i,int j){
   return gsl_matrix_get(_MATRIX,i,j); 
}

//the ith row of the Eigenvector matrix represents the ith Eigenvector
void Matrix::print_row(const char* filename, int ROW){

   FILE * row_file = fopen(filename,"w");

   for(int j=0;j<_N_COLUMN;j++){
      fprintf(row_file,"%3i\t%3.5e\n",(j-((_N_COLUMN-1)/2)),gsl_matrix_get(_MATRIX,ROW,j));
   }
   fclose(row_file);
}


void Matrix::print_to_file(const char* filename){
    
    FILE * matrix_file = fopen(filename,"w");

    for(int i=0;i<_N_ROW;i++){
       for (int j=0;j<_N_COLUMN;j++){
         fprintf(matrix_file,"%3.5e \t", gsl_matrix_get(_MATRIX,i,j));
      }
      fprintf(matrix_file,"\n");
    }
    fclose(matrix_file);
}

void Matrix::print_to_file_list(const char* filename){
    
    FILE * matrix_file = fopen(filename,"w");

    for(int j=0;j<_N_COLUMN;j++){
      for (int i=0;i<_N_ROW;i++){
        fprintf(matrix_file,"%3.5e \n", gsl_matrix_get(_MATRIX,i,j));
      }
    }
    fclose(matrix_file);
}

gsl_matrix Matrix::return_gsl_matrix(){
   return *_MATRIX;
}



