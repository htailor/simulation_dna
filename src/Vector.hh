#ifndef INCLUDE_VECTOR_HH
#define INCLUDE_VECTOR_HH

#include <gsl/gsl_vector.h>

class Vector{
 
public:

 Vector();
 Vector(int n);
 int length();

 double get(int i);

 int Zero_Elements();
 void set(int i, double value);
 void print_to_file(const char* filename);


private:

 int _k;

 gsl_vector *_VECTOR;

};

#endif
