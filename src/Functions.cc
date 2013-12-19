#include "Functions.hh"
#include <sstream>


extern double TOLERANCE;

using namespace std;

double chop (double x){

  double chopped;
  if(abs(x) < TOLERANCE){
    chopped = 0;
  }
  else{
    chopped = x;
  }
  return chopped;
}

double Squared (double x){

  return (x*x); 

}

double Hermite(const int& n, const double& x)
{
    double result;
    int i;
    double a;
    double b;

    // Prepare A and B

    a = 1;
    b = 2*x;

    // Special cases: N=0 or N=1
    if( n==0 )
    {
        result = a;
        return result;
    }
    if( n==1 )
    {
        result = b;
        return result;
    }

    // General case: N>=2
    for(i = 2; i <= n; i++)
    {
        result = 2*x*b-2*(i-1)*a;
        a = b;
        b = result;
    }
    return result;

}

int Factorial(int n){

  int answer = 1, current = 1; 
  while (current <= n) {
    answer *= current;
    current++;
  }
  return answer;
}

double Free_Energy(double x){

  return -log(abs(x));

}

string intToString(int x){

    stringstream ss;
    std::string s;
    ss << x;
    s = ss.str();
    return s;
}


