#ifndef INCLUDE_PARTITION_FUNCTION_HH
#define INCLUDE_PARTITION_FUNCTION_HH

#include <complex>
#include <vector>

using namespace std;

complex<double> EXP_sVt00(int s, int t);
complex<double> EXP_CsVt00(int s, int t);
complex<double> EXP_sVt11(int s, int t);

vector<double> Partition_Function(int nIntact,int nBroken,int Extension);
vector<double> pf_double_frayed(int nBroken1,int nBroken2,int Extension);

double ZNJ_ETA(int Extension,int nj);
double ZNJ_XI(int Extension,int nj);

#endif
