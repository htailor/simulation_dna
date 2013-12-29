#include "Eigenfunctions.hh"
#include "Functions.hh"
#include <cmath>
#include <complex>

using namespace std;

#define PI 3.1415926535897932384626433832795

extern int m;
extern double L;

// Analytical expression of the eigenfunctions of the backbone

complex<double> PSI_S(int s, double eta)
{
	double theta = (s*eta*2*PI)/L;
	double Const = pow(L,-0.5);
	return polar(Const,theta);
}

// Complex function from within the partition function equation

complex<double> g(int s, double eta)
{
	double theta = -(4.0*(double)s*eta*PI)/L;
	return polar(1.0,theta);
}
