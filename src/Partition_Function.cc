#include <cmath>
#include <string>
#include <sstream>
#include <complex>
#include <vector>
#include "Potential.hh"
#include "Matrix.hh"
#include "Vector.hh"
#include "Eigenfunctions.hh"
#include "Eigenvalues.hh"

#define PI 3.1415926535897932384626433832795

using namespace std;

extern double L;
extern double eta_b;
extern int m;
extern double Delta;
extern int _SIZE;
extern int N;

extern int smax;
extern int tmax;
extern int vmax;

extern Vector t00_eval;
extern Matrix t00_evec;

extern Vector t11_eval;
extern Matrix t11_evec;

extern Matrix DoubleIntegral_T10;
extern Matrix DoubleIntegral_T01;

extern vector<vector<complex<double> > > sVt00;
extern vector<vector<complex<double> > > CsVt00;
extern vector<vector<complex<double> > > sVt11;

extern vector<vector<complex<double> > > s1_xi_s2;
extern vector<vector<double> > t_eta_t2;


complex<double> EXP_sVt00(int s, int t)
{
	complex<double> _EXP_sVt00 = (0,0);
	for(int i=-m;i<=m;i++){
		_EXP_sVt00 = _EXP_sVt00 + Delta*PSI_S(s,Delta*i)*exp(-Potential(Delta*i)/2)*t00_evec.get(t,i+m);
	}
	return _EXP_sVt00;
}

complex<double> EXP_CsVt00(int s, int t)
{
	complex<double> _EXP_CsVt00 = (0,0);
	for(int i=-m;i<=m;i++){
		_EXP_CsVt00 = _EXP_CsVt00 + Delta*conj(PSI_S(s,Delta*i))*exp(-Potential(Delta*i)/2)*t00_evec.get(t,i+m);
	}
	return _EXP_CsVt00;
}


complex<double> EXP_sVt11(int s, int t)
{
	complex<double> _EXP_sVt11 = (0,0);
	for(int i=-m;i<=m;i++){
		_EXP_sVt11 = _EXP_sVt11 + Delta*PSI_S(s,Delta*i)*exp(-Potential(Delta*i)/2)*t11_evec.get(t,i+m);
	}
	return _EXP_sVt11;
}

vector<double> Partition_Function(int nIntact,int nBroken,int Extension){
	vector<double> _PF_DATA (2);
	double Sum_PF = 0;
	double Sum_DPF = 0;
	if(nBroken==0){ // Intact Partition Function
		#pragma omp parallel for reduction(+:Sum_PF,Sum_DPF)
		for(int s=-smax;s<=smax;s++){
        	complex<double> s_complex(0.0,(-(4.*PI*s)/L));
			for(int t=0;t<tmax;t++){
				Sum_PF = Sum_PF + pow(lambda_s(s),(N-1))
					*pow(t00_eval.get(t),(N-1))
					*real(
						g(s,Extension*Delta)	
						*sVt00[s+smax][t]
						*sVt00[s+smax][t]
					);
				Sum_DPF = Sum_DPF + pow(lambda_s(s),(N-1))
					*pow(t00_eval.get(t),(N-1))
					*real(
						s_complex
						*g(s,Extension*Delta)
						*pow(sVt00[s+smax][t],2)
					);
			}
		}
		printf("***Intact calculation for extension %3i complete.***\n",Extension);
	}
	else if(nIntact==0){ 
		#pragma omp parallel for reduction(+:Sum_PF,Sum_DPF)
		for(int s=-smax;s<=smax;s++){
        	complex<double> s_complex(0.0,(-(4.*PI*s)/L));
			for(int t=0;t<tmax;t++){
				for(int v=0;v<vmax;v++){
					Sum_PF = Sum_PF + pow(lambda_s(s),(N-1))
						*pow(t00_eval.get(t),(N-nBroken))
						*pow(t11_eval.get(v),(nBroken-1))
						*real(
							g(s,Extension*Delta)
							*sVt00[s+smax][t]
							*sVt11[s+smax][v]
							*DoubleIntegral_T01.get(t,v)
						);
					Sum_DPF = Sum_DPF + pow(lambda_s(s),(N-1))
						*pow(t00_eval.get(t),(N-nBroken))
						*pow(t11_eval.get(v),(nBroken-1))
						*real(
							s_complex
							*g(s,Extension*Delta)
							*sVt00[s+smax][t]
							*sVt11[s+smax][v]
							*DoubleIntegral_T01.get(t,v)
						);
				}
			}
		}
		printf("***Frayed calculation for extension %3i complete.***\n",Extension);
	}
	else{ //Bubble Partition Function
		#pragma omp parallel for reduction(+:Sum_PF,Sum_DPF)
		for(int s=-smax;s<=smax;s++){
        	complex<double> s_complex(0.0,(-(4.*PI*s)/L));
			for(int t=0;t<tmax;t++){
				for(int v=0;v<vmax;v++){
					for(int t2=0;t2<tmax;t2++){
						Sum_PF = Sum_PF + pow(lambda_s(s),(N-1))
						*pow(t00_eval.get(t),(nIntact-1))
						*pow(t00_eval.get(t2),(N-nIntact-nBroken-1))
						*pow(t11_eval.get(v),(nBroken-1))
						*real(
							g(s,Extension*Delta)
							*sVt00[s+smax][t]
							*sVt00[s+smax][t2]
							*DoubleIntegral_T01.get(t2,v)
							*DoubleIntegral_T10.get(v,t)
						);
						Sum_DPF = Sum_DPF + pow(lambda_s(s),(N-1))
						*pow(t00_eval.get(t),(nIntact-1))
						*pow(t00_eval.get(t2),(N-nIntact-nBroken-1))
						*pow(t11_eval.get(v),(nBroken-1))
						*real(
							s_complex
							*g(s,Extension*Delta)
							*sVt00[s+smax][t]
							*sVt00[s+smax][t2]
							*DoubleIntegral_T01.get(t2,v)
							*DoubleIntegral_T10.get(v,t)
						);
					}
				}
			}
		}
		printf("***Bubble calculation for extension %3i complete.***\n",Extension);
	}
	_PF_DATA[0] = Sum_PF;
	_PF_DATA[1] = Sum_DPF;
	return _PF_DATA;
}

vector<double> pf_double_frayed(int nBroken1,int nBroken2,int Extension){
	vector<double> _PF_DATA (2);
	double Sum_PF = 0;
	double Sum_DPF = 0;
	#pragma omp parallel for reduction(+:Sum_PF,Sum_DPF)
	for(int s=-smax;s<=smax;s++){
        complex<double> s_complex(0.0,(-(4.*PI*s)/L));
		for(int v=0;v<vmax;v++){
			for(int v2=0;v2<vmax;v2++){
				for(int t=0;t<tmax;t++){
					Sum_PF = Sum_PF + pow(lambda_s(s),(N-1))
						*pow(t11_eval.get(v),(nBroken1-1))
						*pow(t11_eval.get(v2),(nBroken2-1))
						*pow(t00_eval.get(t),(N-nBroken1-nBroken2-1))
						*real(
							g(s,Extension*Delta)
							*sVt11[s+smax][v]
							*DoubleIntegral_T10.get(v,t)
							*DoubleIntegral_T01.get(t,v2)
							*sVt11[s+smax][v2] // D
						);
					Sum_DPF = Sum_DPF + pow(lambda_s(s),(N-1))
						*pow(t11_eval.get(v),(nBroken1-1))
						*pow(t11_eval.get(v2),(nBroken2-1))
						*pow(t00_eval.get(t),(N-nBroken1-nBroken2-1))
						*real(
							g(s,Extension*Delta)
							*s_complex
							*sVt11[s+smax][v]
							*DoubleIntegral_T10.get(v,t)
							*DoubleIntegral_T01.get(t,v2)
							*sVt11[s+smax][v2] // D
						);
				}
			}
		}
	}
	printf("***Double frayed calculation for extension %3i complete.***\n",Extension);
	_PF_DATA[0] = Sum_PF;
	_PF_DATA[1] = Sum_DPF;
	return _PF_DATA;
}

double ZNJ_ETA(int Extension,int nj){ // This calculates the exact ZNJ. The calculates include an addition nested loop over t'.
	double Sum_ZNJ_ETA = 0;
	#pragma omp parallel for reduction(+:Sum_ZNJ_ETA)
	for(int s1=-smax;s1<=smax;s1++){
		for(int t1=0;t1<tmax;t1++){
			for(int t2=0;t2<tmax;t2++){
				Sum_ZNJ_ETA = Sum_ZNJ_ETA + pow(lambda_s(s1),(N-1))
                    *pow(t00_eval.get(t2),(N-nj))
                    *pow(t00_eval.get(t1),(nj-1))
                    *real(
                          g(s1,Extension*Delta)
                          *sVt00[s1+smax][t1]
                          *sVt00[s1+smax][t2]
                          *t_eta_t2[t1][t2]
                          );
			}
		}
	}
	return Sum_ZNJ_ETA;
}

double ZNJ_XI(int Extension,int nj){ // This calculates the exact ZNJ. The calculates include an addition nested loop over t'.
	double Sum_ZNJ_XI = 0;
	#pragma omp parallel for reduction(+:Sum_ZNJ_XI)
	for(int s1=-smax;s1<=smax;s1++){
		for(int s2=-smax;s2<=smax;s2++){
			for(int t1=0;t1<tmax;t1++){
				Sum_ZNJ_XI = Sum_ZNJ_XI + pow(t00_eval.get(t1),(N-1))
				*pow(lambda_s(s1),(N-nj))
				*pow(lambda_s(s2),(nj-1))
				*real(  
					g(s1,Extension*Delta)
					*sVt00[s1+smax][t1]
					*CsVt00[s2+smax][t1]
					*s1_xi_s2[s2+smax][s1+smax]
				);
			}
		}
	}
	return Sum_ZNJ_XI;
}
