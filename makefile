#///////////////////////////////////////////////////
#//
#//  Program created by Paul Murphy
#//  
#//  This is a C program version of the 
#//  the sample SCF Calculation found in Appendix B
#//  of Modern Quantum Chemistry, Introduction to 
#//  Advanced Electronic Structure Theory by 
#//  A. Szabo and N. Ostlund (originally written in 
#//  Fortran IV). Comments and feedback are welcome
#//  by email - juansanshoo@hotmail.co.uk.
#//
#//  MAKEFILE
#//  --------
#//  Version History: 1.0 - 31/12/2015
#//
#///////////////////////////////////////////////////
#DEST = ./
#OBJS = main.o calc_integrals.o format_integrals.o perform_scf.o
scf_basic_version_1_0: main.c calc_integrals.c format_integrals.c perform_scf.c
	gcc -o scf_basic_version_1_0 main.c calc_integrals.c format_integrals.c perform_scf.c -I. -lm
