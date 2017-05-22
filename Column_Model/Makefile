#
# This is the makefile .
# To compile the code, just type 'make' and then type './Run_1DCol'
#

#COMPILER = gfortran
COMPILER = ifort
#COMPILER = g95
COMPOPTS = -r8 -g
#COMPOPTS = -g -ffixed-line-length-none

column_lkPro.o: column_lkPro.f90
	$(COMPILER)  column_lkPro.f90 -o Run_1Dcol
#	$(COMPILER) column_lkPro.f90 -o -g -traceback Run_1Dcol 

# To clean previous outputs for a new run, it is recommended that you type 'make clean' first
# NOTE: This option removes all output files 
clean:
	rm -f *.o fort.* *~ *.mod *.exe *.out Main core


