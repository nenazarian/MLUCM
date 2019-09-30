#     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#     +  Created and modified by N. Nazarian, Singapore-MIT/UNSW       +
#     +  Jan 2017-Sept 2019                                            +
#     +                        email: n.nazarian@unsw.edu.au/		   +
#	  +								  nazarian@mit.edu                 +
#     +                             /nenazarian@gmail.com              +
#     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#     Goal: This is the makefile to compile and run the 1D column model. 
#	  Procedure: To compile the code, type 'make' which will create the+
#                executable file "Run_1DCol". 						   +
#				 Afterward, type './Run_1DCol' to run the column model.+
#     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# The compiler used here is ifort. This can be changed based on the working 
# FORTRAN compiler
#COMPILER = gfortran
COMPILER = ifort
#COMPILER = g95
COMPOPTS = -r8 -g

column.o: column.f90
		$(COMPILER)  column.f90 -o Run_1Dcol

# To clean previous outputs for a new compilation, it is recommended   +
# that you type 'make clean' first. 								   +
# !!!! NOTE: This option removes all output files !!!				   +
clean:
	rm -f Run_1Dcol *.o fort.* *~ *.mod *.exe *.out Main core
#     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
