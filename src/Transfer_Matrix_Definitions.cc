#include "Transfer_Matrix_Definitions.hh"
#include "Functions.hh"
#include <cmath>

extern double eta_b;

using namespace std;

double g0(double eta)
{
   double _g0;
   if(abs(eta) <= eta_b){
      _g0 = 1;
   }
   else{
      _g0 = 0;
   }
   return _g0;
}

double g1(double eta)
{
   double _g1;
   if(abs(eta) <= eta_b){
      _g1 = 0;
   }
   else{
      _g1 = 1;
   }
   return _g1;
}

double T(double ea, double eb)
{
   return exp(-Squared(ea-eb));
}

double T_hat(double na, double nb)
{
   return exp(-Squared(na-nb))*exp(-(Potential(na) + Potential(nb))/2); 
}

double T_hat00(double na, double nb)
{
   return g0(na)*g0(nb)*T_hat(na,nb);
}

double T_hat01(double na, double nb)
{
   return g0(na)*g1(nb)*T_hat(na,nb);
}

double T_hat10(double na, double nb)
{
   return g1(na)*g0(nb)*T_hat(na,nb);
}

double T_hat11(double na, double nb)
{
   return g1(na)*g1(nb)*T_hat(na,nb);
}

