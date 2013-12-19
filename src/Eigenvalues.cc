#include <cmath>
#include "Eigenvalues.hh"
#include <stdlib.h>

#define PI 3.1415926535897932384626433832795

extern double L;

double lambda_s(int s)
{
   return pow(PI,0.5)*exp(-pow((PI*s)/L,2.));
}
