#include "Vector.hh"
#include <cmath>

using namespace std;

Vector::Vector(){}

Vector::Vector(int n){
   
   _k = n;
   _VECTOR = gsl_vector_alloc(_k);
}

int Vector::length(){
   return _k;
}

void Vector::set(int i, double value){
   gsl_vector_set(_VECTOR,i,value); 
}

double Vector::get(int i){
   return gsl_vector_get(_VECTOR,i); 
}

void Vector::print_to_file(const char* filename){
    
    FILE * vector_file = fopen(filename,"w");

    for(int i=0;i<_k;i++){
      fprintf(vector_file,"%3.5e \n", gsl_vector_get(_VECTOR,i));
    }
    fclose(vector_file);

}


int Vector::Zero_Elements(){

    int count=0;
    for(int i=0 ; i<_k ; i++){
        if( abs(gsl_vector_get(_VECTOR,i)) == 0 ){ 
            count = count + 1; 
        }
    } 
    return count;
}
