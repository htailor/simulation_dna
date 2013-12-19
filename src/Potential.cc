#include "Potential.hh"
#include "Functions.hh"
#include <cmath>
#include <iostream>

using namespace std;

extern double Beta;
extern double e0;
extern double eta_d;

double Potential(double _n){
    
  double PotentialValue = -((Beta*e0)/pow(1 + 2.*(Squared(_n)/Squared(eta_d)),3.));
  return PotentialValue;
}

