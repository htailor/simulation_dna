#include "Matrix_Functions.hh"
#include "Functions.hh"
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <cmath>

using namespace std;

Matrix Eigenvectors(Matrix M){

    int _k = M.n_row();

    gsl_matrix *_M = gsl_matrix_alloc(_k,_k);

    for(int i=0;i<_k;i++){
       for (int j=0;j<_k;j++){
         gsl_matrix_set(_M,i,j,M.get(i,j));
      }
    }

    gsl_matrix *EVEC = gsl_matrix_alloc(_k,_k);
    gsl_vector *DUMMY = gsl_vector_alloc(_k);

    gsl_eigen_symmv_workspace *TEMP_EVEC = gsl_eigen_symmv_alloc(_k);
    gsl_eigen_symmv(_M,DUMMY,EVEC,TEMP_EVEC);

    gsl_eigen_symmv_sort(DUMMY,EVEC,GSL_EIGEN_SORT_ABS_DESC);
    gsl_eigen_symmv_free(TEMP_EVEC);

    gsl_matrix_transpose(EVEC);
   
    Matrix eigenvectors = Matrix(_k,_k);

    for(int i=0;i<_k;i++){
       for (int j=0;j<_k;j++){
           eigenvectors.set(i,j,gsl_matrix_get(EVEC,i,j));
           //eigenvectors.set(i,j,chop(gsl_matrix_get(EVEC,i,j)));
       }
    }
    return eigenvectors;
}

Vector Eigenvalues(Matrix M){

    int _k = M.n_row();

    gsl_matrix *_M = gsl_matrix_alloc(_k,_k);

    for(int i=0;i<_k;i++){
       for (int j=0;j<_k;j++){
         gsl_matrix_set(_M,i,j,M.get(i,j));
      }
    }

    gsl_matrix *DUMMY = gsl_matrix_alloc(_k,_k);
    gsl_vector *EVAL = gsl_vector_alloc(_k);

    gsl_eigen_symmv_workspace *_TEMP_ = gsl_eigen_symmv_alloc(_k);
    gsl_eigen_symmv(_M,EVAL,DUMMY,_TEMP_);
    
    gsl_eigen_symmv_sort(EVAL,DUMMY,GSL_EIGEN_SORT_ABS_DESC);

    gsl_eigen_symmv_free(_TEMP_);

    Vector eigenvalues = Vector(_k);

    for(int i=0;i<_k;i++){
         eigenvalues.set(i,chop(gsl_vector_get(EVAL,i)));
    }
    return eigenvalues;

}
