#include "Matrix.hh"
#include "Vector.hh"
#include "Potential.hh"
#include "Matrix_Functions.hh"
#include "Partition_Function.hh"
#include "Transfer_Matrix_Definitions.hh"
#include "Functions.hh"
#include "Eigenfunctions.hh"
#include "Eigenvalues.hh"
#include "Menu.hh"
#include <cmath>
#include <getopt.h>
#include <ctime>
#include <complex>
#include <vector>
#include <stdio.h>
#include <omp.h>
#include <iostream>
#include <sstream>

using namespace std;

// Global Parameters//

double Beta = 38.6817; // eV^-1
double eta_d = 20.; // diameter of dna AA
double e0 = 0;

int umax=0;
int umin;
int N=0; 
int m=0;
int MATRIX_SIZE;

int nIntact;
int nBroken;

int smax=0;
int tmax=0;
int pmax=0;
int vmax=0;

double L=0;
double eta_b=0;
double kappa=0;
double sigma=0;
double kappa_sigma_r=0;
double TOLERANCE;
double Delta;

Matrix t;
Matrix t_evec;
Vector t_eval;

Matrix t00;
Matrix t00_evec;
Vector t00_eval;

Matrix t11;
Matrix t11_evec;
Vector t11_eval;

Matrix DoubleIntegral_T10;
Matrix DoubleIntegral_T01;

vector<vector<complex<double> > > sVt00;
vector<vector<complex<double> > > CsVt00;
vector<vector<complex<double> > > sVt11;
vector<vector<double> > t_eta_t2;
vector<vector<double> > t11_t00;
vector<vector<complex<double> > > s1_xi_s2;


void NucleationRun(int intact, int broken){

	/*
		Calculates the partition function for a nucleated state specified by the number of intact and broken base-pairs.
	*/

	// Create files for free energy and partition function
	string FE_filename  = "FreeEnergy_" + intToString(intact) + "_" + intToString(broken) + ".out";
	string PF_filename  = "PartitionFunction_" + intToString(intact) + "_" + intToString(broken) + ".out";

	Matrix FE = Matrix(umax-umin+1,2);
	Matrix PF = Matrix(umax-umin+1,2);
	
	// Create files for the force (first derivative of the free energy and partition function)
	string dFE_filename = "dFreeEnergy_" + intToString(intact) + "_" + intToString(broken) + ".out";
	string dPF_filename = "dPartitionFunction_" + intToString(intact) + "_" + intToString(broken) + ".out";
  
	Matrix dFE = Matrix(umax-umin+1,2);
	Matrix dPF = Matrix(umax-umin+1,2);

	// Create files for the mean axial displacements
	string ZNJ_xi_filename = "znj_xi.out";
	string ZNJ_eta_filename = "znj_eta.out";

	Matrix ZNJ_xi = Matrix(umax-umin+1,N);
	Matrix ZNJ_eta = Matrix(umax-umin+1,N);

	string MXD_xi_filename = "MeanAxialDisp_xi.out";
	string MXD_eta_filename = "MeanAxialDisp_eta.out";

	Matrix MXD_xi = Matrix(umax-umin+1,N);
	Matrix MXD_eta = Matrix(umax-umin+1,N);

	string Average_X_filename = "Average_X.out";
	string Average_Y_filename = "Average_Y.out";

	Matrix AverageX = Matrix(umax-umin+1,N);
	Matrix AverageY = Matrix(umax-umin+1,N);
	
	for(int i=umin;i<=umax;i++){

		string ext_status_create = "touch RUNNING_EXT_" + intToString(i) + ".status";
		string ext_status_destroy = "rm -f RUNNING_EXT_" + intToString(i) + ".status";
		system(ext_status_create.c_str());

		vector<double> PARTITION_FUNCTION_DATA(Partition_Function(intact,broken,i));
		double PARTITION_FUNCTION_VALUE = PARTITION_FUNCTION_DATA[0];
		double dPARTITION_FUNCTION_VALUE = PARTITION_FUNCTION_DATA[1];
		double dFreeEnergy = -(1/PARTITION_FUNCTION_VALUE)*dPARTITION_FUNCTION_VALUE;

		PF.set(i-umin,1,PARTITION_FUNCTION_VALUE);
		PF.set(i-umin,0,Delta*i);

		FE.set(i-umin,1,Free_Energy(PARTITION_FUNCTION_VALUE));
		FE.set(i-umin,0,Delta*i);
			
		dPF.set(i-umin,1,dPARTITION_FUNCTION_VALUE);
		dPF.set(i-umin,0,Delta*i);	

		dFE.set(i-umin,1,dFreeEnergy);
		dFE.set(i-umin,0,Delta*i);

		// Print the data to files
    	PF.print_to_file(PF_filename.c_str());
    	dPF.print_to_file(dPF_filename.c_str());

    	FE.print_to_file(FE_filename.c_str());
		dFE.print_to_file(dFE_filename.c_str());

		/*if(broken==0){  // Intact State Only - Calculate MXD
			for(int j=1;j<=N;j++){ // j represents the base-pair, i refers the extension applied.
				string mxd_status_create = "touch RUNNING_MXD_" + intToString(j) + ".status";
				string mxd_status_destroy = "rm -f RUNNING_MXD_" + intToString(j) + ".status";
				system(mxd_status_create.c_str());

        		double ZNJ_XI_VALUE = ZNJ_XI(i,j);
				double ZNJ_ETA_VALUE = ZNJ_ETA(i,j); // refer to BasePairCalculations.cc for definition.

      			double MXD_XI_VALUE = ZNJ_XI_VALUE/PARTITION_FUNCTION_VALUE;
	        	double MXD_ETA_VALUE = ZNJ_ETA_VALUE/PARTITION_FUNCTION_VALUE;

        		double AVERAGE_X_VALUE = 0.5*(MXD_XI_VALUE+MXD_ETA_VALUE);
        		double AVERAGE_Y_VALUE = 0.5*(MXD_XI_VALUE-MXD_ETA_VALUE);

	        	ZNJ_xi.set(i,j-1,ZNJ_XI_VALUE);
        		MXD_xi.set(i,j-1,ZNJ_XI_VALUE/PARTITION_FUNCTION_VALUE);

				ZNJ_eta.set(i,j-1,ZNJ_ETA_VALUE);
				MXD_eta.set(i,j-1,ZNJ_ETA_VALUE/PARTITION_FUNCTION_VALUE);

	        	AverageX.set(i,j-1,AVERAGE_X_VALUE);
        		AverageY.set(i,j-1,AVERAGE_Y_VALUE);

				system(mxd_status_destroy.c_str());
			}
		
		}*/
		system(ext_status_destroy.c_str());    

	}

	if(broken==0){

		//ZNJ_xi.print_to_file(ZNJ_xi_filename.c_str());
		//ZNJ_eta.print_to_file(ZNJ_eta_filename.c_str());

		//MXD_eta.print_to_file(MXD_eta_filename.c_str());
		//MXD_xi.print_to_file(MXD_xi_filename.c_str());

		//AverageX.print_to_file(Average_X_filename.c_str());
		//AverageY.print_to_file(Average_Y_filename.c_str());

		printf("\n***INTACT STATE COMPLETE***\n\n");
		system("mv *.out ./results/Intact");
	}
	else if(intact==0){
		printf("\n***FRAYED STATE COMPLETE %2i/%2i***\n\n",intact,broken);
		system("mv *.out ./results/Frayed");
	}
	else{
		printf("\n***BUBBLE STATE COMPLETE %2i/%2i***\n\n",intact,broken);
		system("mv *.out ./results/Bubble");
	}	

}

// Methods for double frayed states from the ends of the DNA structure
void NucleationDoubleFrayed(int broken1, int broken2){

	string df_status_create = "touch RUNNING_DOUBLE_FRAYED_" + intToString(broken1) + "_" + intToString(broken2) + ".status";
	string df_status_destroy = "rm -f RUNNING_DOUBLE_FRAYED_" + intToString(broken1) + "_" + intToString(broken2) + ".status";
	system(df_status_create.c_str());

	string FE_filename  = "DoubleFrayed_FreeEnergy_" + intToString(broken1) + "_" + intToString(broken2) + ".out";
	string PF_filename  = "DoubleFrayed_PartitionFunction_" + intToString(broken1) + "_" + intToString(broken2) + ".out";

	Matrix FE = Matrix(umax-umin+1,2);
	Matrix PF = Matrix(umax-umin+1,2);

	string dFE_filename = "DoubleFrayed_dFreeEnergy_" + intToString(broken1) + "_" + intToString(broken2) + ".out";
	string dPF_filename = "DoubleFrayed_dPartitionFunction_" + intToString(broken1) + "_" + intToString(broken2) + ".out";
  
	Matrix dFE = Matrix(umax-umin+1,2);
	Matrix dPF = Matrix(umax-umin+1,2);

	for(int i=umin;i<=umax;i++){

		string ext_status_create = "touch RUNNING_DF_EXT_" + intToString(i) + ".status";
		string ext_status_destroy = "rm -f RUNNING_DF_EXT_" + intToString(i) + ".status";
		system(ext_status_create.c_str());

		vector<double> PARTITION_FUNCTION_DATA(pf_double_frayed(broken1,broken2,i));
		double PARTITION_FUNCTION_VALUE = PARTITION_FUNCTION_DATA[0];
		double dPARTITION_FUNCTION_VALUE = PARTITION_FUNCTION_DATA[1];
		double dFreeEnergy = -(1/PARTITION_FUNCTION_VALUE)*dPARTITION_FUNCTION_VALUE;

		PF.set(i-umin,1,PARTITION_FUNCTION_VALUE);
		PF.set(i-umin,0,Delta*i);

		FE.set(i-umin,1,Free_Energy(PARTITION_FUNCTION_VALUE));
		FE.set(i-umin,0,Delta*i);
			
		dPF.set(i-umin,1,dPARTITION_FUNCTION_VALUE);
		dPF.set(i-umin,0,Delta*i);	

		dFE.set(i-umin,1,dFreeEnergy);
		dFE.set(i-umin,0,Delta*i);

        PF.print_to_file(PF_filename.c_str());
    	dPF.print_to_file(dPF_filename.c_str());

    	FE.print_to_file(FE_filename.c_str());
		dFE.print_to_file(dFE_filename.c_str());

		system("mv *.out ./results/Frayed");

		system(ext_status_destroy.c_str());    
	}
	system(df_status_destroy.c_str());
	printf("\n***DOUBLE FRAYED STATE COMPLETE %2i/%2i***\n\n",broken1,broken2);
}

int main(int argc,char *argv[]) 
{
   //omp_set_num_threads(8);
   Menu(argc,argv);    

   system("rm -f *.status");
   system("touch SIM_STARTED.status");
   system("touch RUNNING.status");

   time_t NUCLEATION_START,NUCLEATION_FINISH;	// Start timer to calculate total runtime
   time (&NUCLEATION_START);
   printf("*******************************************************************************************\n");
   printf("Starting Nucleation Simulation.....\n");
   printf("*******************************************************************************************\n\n");

   printf("Creating data for base-pair potential...\n");

   Matrix PotentialData = Matrix(MATRIX_SIZE,2);

   for(int i=-m;i<=m;i++){
   	  double POTENTIAL_VALUE = Potential(i*Delta);
   	  PotentialData.set(i+m,1,POTENTIAL_VALUE);
   	  PotentialData.set(i+m,0,Delta*i);
   }

   PotentialData.print_to_file("Potential_Data.out");
   system("mv Potential_Data.out ./results");
   printf("Complete.\n\n");

///////////////////////////////////////////////////////////////////////////////////////

   ///STAGE 1 - CREATE TRANSFER MATRICES T,T00,T11///

///////////////////////////////////////////////////////////////////////////////////////

   /*
   Create the transfer matrices T,T00,T11 using the Transfer Matrix Functions (These functions
   contain the mathematical definitions and hence calculate the values for each element in the
   matrix
   */

   printf("Initialising Transfer Matrices T, T00, T11.....");

   Matrix t00 (MATRIX_SIZE,MATRIX_SIZE);    
   Matrix t11 (MATRIX_SIZE,MATRIX_SIZE);   

   for(int i=0;i<=2*m;i++){
      for(int j=0;j<=2*m;j++){
         t00.set(i,j,Delta*T_hat00(Delta*(i-m),Delta*(j-m)));  // T_hat00 definition is in Transfer_Matrix_Definitions.cc   
         t11.set(i,j,Delta*T_hat11(Delta*(i-m),Delta*(j-m)));  // T_hat11 definition is in Transfer_Matrix_Definitions.cc  
      }
   }

   
   printf("COMPLETE.\n\n");  
   

///////////////////////////////////////////////////////////////////////////////////////
 
   ///STAGE 2 - CALCULATE EIGENVALUES AND EIGENVECTORS OF MATRICES T,T00,T11///

///////////////////////////////////////////////////////////////////////////////////////

   /*
   STAGE 2 solves the eigensystem for matrices T,T00 and T11. The Eigenvalues are stored in
   a vector and the eigenvectors are stored in a matrix. 
   */

   printf("Calculating Eigenvalues and Eigenvectors of T, T00 and T11.....");  

   time_t Eigenvalue_Problem_START,Eigenvalue_Problem_FINISH;
   time (&Eigenvalue_Problem_START);

   t00_eval = Eigenvalues(t00);
   t00_evec = Eigenvectors(t00);

   t11_eval = Eigenvalues(t11);
   t11_evec = Eigenvectors(t11);

   for(int i=0;i<MATRIX_SIZE;i++){
      for(int j=0;j<MATRIX_SIZE;j++){
         double NewElement_T00 = pow(Delta,-0.5)*t00_evec.get(i,j);
         t00_evec.set(i,j,NewElement_T00);
         double NewElement_T11 = pow(Delta,-0.5)*t11_evec.get(i,j);
         t11_evec.set(i,j,NewElement_T11);
      }
   }

   smax = m;
   tmax = MATRIX_SIZE-t00_eval.Zero_Elements();
   vmax = MATRIX_SIZE-t11_eval.Zero_Elements();
   
   printf("\n\n------------------------------------------------\n");
   printf("smax: %4i tmax: %4i vmax: %4i",smax,tmax,vmax);
   printf("\n------------------------------------------------\n\n");

   time (&Eigenvalue_Problem_FINISH);
   printf("COMPLETED in %6.2f seconds.\n\n", difftime(Eigenvalue_Problem_FINISH,Eigenvalue_Problem_START));

///////////////////////////////////////////////////////////////////////////////////////

   ///STAGE 3 - CALCULATE DOUBLE INTEGRALS FOR T10 and T01 Matrices///

///////////////////////////////////////////////////////////////////////////////////////

   // Calculates all the integral sums for the partition function calculations

   sVt00 = vector<vector<complex<double> > > (2*smax+1, vector<complex<double> > (tmax) );
   CsVt00 = vector<vector<complex<double> > > (2*smax+1, vector<complex<double> > (tmax) );
   sVt11 = vector<vector<complex<double> > > (2*smax+1, vector<complex<double> > (vmax) );

   s1_xi_s2 = vector<vector<complex<double> > > (2*smax+1, vector<complex<double> > (2*smax+1) );
   t_eta_t2 = vector<vector<double > > (tmax, vector<double> (tmax) );

   DoubleIntegral_T10 = Matrix(vmax,tmax);
   DoubleIntegral_T01 = Matrix(tmax,vmax);

   time_t DI_START, DI_FINISH;
   time (&DI_START);

   printf("Pre-calculating Standard Integrals.....");

/// Pre Calculating the Double Intergrals ///

   #pragma omp parallel for
   for(int t=0;t<tmax;t++){
      for(int v=0;v<vmax;v++){
    	 double _DI_T10 = 0;
         double _DI_T01 = 0;
         for(int i=-m;i<=m;i++){
            for(int j=-m;j<=m;j++){
                _DI_T10 = _DI_T10 + Delta*Delta*t11_evec.get(v,i+m)*T_hat10(Delta*i,Delta*j)*t00_evec.get(t,j+m);
                _DI_T01 = _DI_T01 + Delta*Delta*t00_evec.get(t,i+m)*T_hat01(Delta*i,Delta*j)*t11_evec.get(v,j+m);
            }
         }
         DoubleIntegral_T10.set(v,t,_DI_T10);
         DoubleIntegral_T01.set(t,v,_DI_T01);
      }
   }

   #pragma omp parallel for
   for(int t=0;t<tmax;t++){
      for(int t2=0;t2<tmax;t2++){
    	 double Integral_1=0;
    	 for(int i=-m;i<=m;i++){
    	    Integral_1 = Integral_1 + Delta*t00_evec.get(t,i+m)*Delta*i*t00_evec.get(t2,i+m);
    	 }
    	 t_eta_t2[t][t2] = Integral_1;
      }
   }

   #pragma omp parallel for
   for(int s1=-smax;s1<=smax;s1++){
      for(int s2=-smax;s2<=smax;s2++){
    	 complex<double> Integral_2 = (0,0);
    	 for(int i=-m;i<=m;i++){
    	    Integral_2 = Integral_2 + Delta*PSI_S(s1,i*Delta)*(Delta*i)*PSI_S(s2,i*Delta);
    	 }
    	 s1_xi_s2[s1+smax][s2+smax] = Integral_2;
      }
   }

   printf("COMPLETE.\n\nPre-calculating Complex Integrals.....COMPLETE\n\n");

   #pragma omp parallel for
   for(int s=-smax;s<=smax;s++){
      for(int t=0;t<tmax;t++){
          sVt00[s+smax][t] = EXP_sVt00(s,t);
          CsVt00[s+smax][t] = EXP_CsVt00(s,t);
      }
   }

   #pragma omp parallel for
   for(int s=-smax;s<=smax;s++){
      for(int v=0;v<vmax;v++){
     	 sVt11[s+smax][v] = EXP_sVt11(s,v);
      }
   }

   time (&DI_FINISH);
   printf("COMPLETED in %6.2f seconds.\n\n",difftime(DI_FINISH,DI_START));

///////////////////////////////////////////////////////////////////////////////////////
 
   ///STAGE 4 - CALCULATE PARTITION FUNCTION///

///////////////////////////////////////////////////////////////////////////////////////

   printf("Calculating Partition Functions.....\n\n");   

   time_t PF_TOTAL_START, PF_TOTAL_FINISH;
   time (&PF_TOTAL_START);

   //This calculates the Intact State
   system("touch RUNNING_INTACT_STATE.status");
   NucleationRun(0,0);
   system("rm -f RUNNING_INTACT_STATE.status");
   system("touch INTACT_STATE_COMPLETED.status");

   //This calculates all Frayed States
   system("touch RUNNING_FRAYED_STATE.status");
   for(nBroken=1;nBroken<=3;nBroken++){
      NucleationRun(0,nBroken);
   }
   system("rm -f RUNNING_FRAYED_STATE.status");
   system("touch FRAYED_STATE_COMPLETED.status");

   //Calculates the bubbling of three base-pairs from the end of DNA

   system("touch RUNNING_BUBBLE_STATE.status");

   if(N>=9){
       nIntact = 1;
       for(nBroken=1;nBroken<=3;nBroken++){
           NucleationRun(nIntact,nBroken);
       }
   }

   system("rm -f RUNNING_BUBBLE_STATE.status");
   system("touch BUBBLE_STATE_COMPLETED.status");


 //Calculates Double Frayed States, fraying at both ends of DNA
   system("touch RUNNING_DOUBLE_FRAYED_STATE.status");

   int nb = 0;	
   if(N%2==0){ // N is even
	nb = (N/2)-1;
   }
   else{ // N is odd
        nb = (N-1)/2;
   }
   
   if(N>7){nb=3;} // Sets maximum fraying to 3 base-pairs at each end

   int counter = 0;
   for(int i=1;i<=nb;i++){
       for(int j=1+counter;j<=nb;j++){
           NucleationDoubleFrayed(i,j);
       }
       counter = counter +1;
   }

   system("rm -f RUNNING_DOUBLE_FRAYED_STATE.status");
   system("touch DOUBLE_FRAYED_STATE_COMPLETED.status");

   time(&PF_TOTAL_FINISH);
   printf("COMPLETED in %6.2f seconds\n\n",difftime(PF_TOTAL_FINISH,PF_TOTAL_START));


///////////////////////////////////////////////////////////////////////////////////////
 
   ///STAGE 5 - CREATE OUTPUT FILES///

///////////////////////////////////////////////////////////////////////////////////////

  
   //This section writes all the data and logs to files. The data output includes
   //eigenvalue, eigenvectors, partition functions, free energies, Transfer matrices.
   
   //cout << "Writing output data to files....\n";
   //OutputDataToFiles();

   time (&NUCLEATION_FINISH);
   printf("---------------------------\n");
   printf("NUCLEATION RUNTIME: %6.2f\n",difftime(NUCLEATION_FINISH,NUCLEATION_START));
   printf("---------------------------\n\n");

   system("rm -f RUNNING.status");
   system("touch SIM_FINISHED.status");
    
   return 0;
}
