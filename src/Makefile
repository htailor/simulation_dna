all: Nucleation

INCLUDE=/usr/include/gsl
LIBS=-lgsl -lgslcblas -lm

COMPILE=g++ -fopenmp -c -g -I$(INCLUDE)
LINK=g++ -fopenmp

all: Nucleation

Nucleation : Nucleation.o Potential.o Matrix.o Vector.o Matrix_Functions.o Eigenvalues.o Eigenfunctions.o Transfer_Matrix_Definitions.o Partition_Function.o Functions.o Menu.o
	echo "<**Linking**> Nucleation"
	$(LINK) Nucleation.o Potential.o Matrix.o Vector.o Matrix_Functions.o Eigenvalues.o Eigenfunctions.o Transfer_Matrix_Definitions.o Partition_Function.o Functions.o Menu.o $(LIBS) -o Nucleation

Nucleation.o : Nucleation.cc 
	echo "<**Compiling**> Nucleation.cc"
	$(COMPILE) Nucleation.cc

Potential.o : Potential.cc Potential.hh
	echo "<**Compiling**> Potential.cc"
	$(COMPILE) Potential.cc

Transfer_Matrix_Definitions.o : Transfer_Matrix_Definitions.cc Transfer_Matrix_Definitions.hh
	echo "<**Compiling**> Transfer_Matrix_Definitions.cc"
	$(COMPILE) Transfer_Matrix_Definitions.cc

Matrix.o : Matrix.cc Matrix.hh
	echo "<**Compiling**> Matrix.cc"
	$(COMPILE) Matrix.cc

Vector.o : Vector.cc Vector.hh
	echo "<**Compiling**> Vector.cc"
	$(COMPILE) Vector.cc

Matrix_Functions.o : Matrix_Functions.cc Matrix_Functions.hh
	echo "<**Compiling**> Matrix_Functions.cc"
	$(COMPILE) Matrix_Functions.cc

Partition_Function.o : Partition_Function.cc Partition_Function.hh
	echo "<**Compiling**> Partition_Function.cc"
	$(COMPILE) Partition_Function.cc

Functions.o : Functions.cc Functions.hh
	echo "<**Compiling**> Functions.cc"
	$(COMPILE) Functions.cc

Eigenfunctions.o : Eigenfunctions.cc Eigenfunctions.hh
	echo "<**Compiling**> Eigenfunctions.cc"
	$(COMPILE) Eigenfunctions.cc

Eigenvalues.o : Eigenvalues.cc Eigenvalues.hh
	echo "<**Compiling**> Eigenvalues.cc"
	$(COMPILE) Eigenvalues.cc

Menu.o : Menu.cc Menu.hh
	echo "<**Compiling**> Menu.cc"
	$(COMPILE) Menu.cc
