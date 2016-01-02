///////////////////////////////////////////////////
//
//	Program created by Paul Murphy
//	
//	This is a C program version of the 
//	the sample SCF Calculation found in Appendix B
//	of Modern Quantum Chemistry, Introduction to 
//	Advanced Electronic Structure Theory by 
//	A. Szabo and N. Ostlund (originally written in 
//	Fortran IV). Comments and feedback are welcome
//	by email - juansanshoo@hotmail.co.uk.
//
//	MAIN.C
//	------
//	Version History: 1.0 - 31/12/2015
//
///////////////////////////////////////////////////



#include "stdio.h"
#include "scf_globals.h"

//0 = none; 1 = only converged results; 2 = every iteration
PRINT_LEVELS print_level = PRINT_ALL;

//Select the n in STO-nG basis set. Hard code for now.
UINT32 sto_ng = 3;

double bond_len = 1.4632; //bond length
double bond_len_sq; //bond length squared
double zetaA = 2.0925; //slater orbital exponent for atom A
double zetaB = 1.24; //slater orbital exponent for atom B
double za = 2.0; //nuclear charge of atom A
double zb = 1.0; //nuclear charge of atom B



void main(void)
{

	printf("Welcome to the C version of the Szabo/Ostlund SCF software from the book Modern Quantum Chemistry.\n");
	printf("\n*****Version 1.0. C version created by Paul Murphy pm133@hotmail.co.uk*****\n");

	if(print_level == PRINT_ALL)
	{
		printf("\n\nPrint option selected: ALL\n");
	}
	else
	{
		printf("\n\nPrint option selected %d\n", print_level);
	}
	
    printf("\nSelected Basis Set is STO-%dG for nuclei %.1f and %.1f\n", sto_ng, za, zb);
	
	//Calculate the one and two electron integrals
	calc_integrals();

	//Put all integrals in nice arrays simply for illustration
	format_integrals();

	//Finally perform the SCF calculation
	perform_scf();

	printf("\n\n**********End of Calculation**********\n");
}
