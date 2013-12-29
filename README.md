DNA Simulation
==============

Partition function calculations for DNA using the Transfer Matrix Method.

Requirements
============

* Linux OS
* Mathematica (data analysis)

Compiling Requirements
======================

* GCC 4.7
* OpenMP
* GSL Scientific Library

Compiling
=========

Run the script `COMPILE.sh` to compile. This will create a new binary file called `Nucleation` with the older version being renamed to `Nucleation.bak`. The source code and makefile are located in the ./src directory. 

Instructions
============

The table highlights the parameters needed to run the calculations

| Parameter | Cmd Flag |              Type              |                          Notes                         |
|:---------:|:--------:|:------------------------------:|:------------------------------------------------------:|
| L         | --L      | &lt;int&gt;  or &lt;double&gt; | Specifies the range of dimensionless space            |
| m         | --m      |           &lt;int&gt;          | Sets the delta                                         |
| N         | --N      |           &lt;int&gt;          | Sets the number of base-pairs in the structure         |
| eta_b     | --eta_b  |         &lt;double&gt;         | Sets the elastic range in the potential                |
| kappa     | --kappa  |         &lt;double&gt;         | Sets the spring constant of the backbones               |
| sigma     | --sigma  |         &lt;double&gt;         | Sets the spring constant of the base-pairs             |
| umin      | --umin   |           &lt;int&gt;          | Starting position of the extension (number of delta's) |
| umax      | --umax   |           &lt;int&gt;          | End position of the extension (number of delta's)      |


There are two methods to run the calculations. The first method involves running the binary file from the command line

```
$ ./Nucleation
```

This will initiate a user interface asking for parameters needed for the calculations. The second method bypasses the user interface and allows the passing of parameters from the command line 

```
$ ./Nucleation --L <INT> --m <INT> --N <INT> --eta_b <DOUBLE> --kappa <DOUBLE> --sigma <DOUBLE> --umin <INT> --umax <INT>
```

Results from the calculations are written to the directory `./results` categorised into sub-directories for each state (Intact, Frayed and Bubble). After the results are generated run the mathematica books in the directories to perform the analysis and plots.
