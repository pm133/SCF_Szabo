///////////////////////////////////////////////////
//
//  Program created by Paul Murphy
//  
//  This is a C program version of the 
//  the sample SCF Calculation found in Appendix B
//  of Modern Quantum Chemistry, Introduction to 
//  Advanced Electronic Structure Theory by 
//  A. Szabo and N. Ostlund (originally written in 
//  Fortran IV). Comments and feedback are welcome
//  by email - juansanshoo@hotmail.co.uk.
//
//  FORMAT_INTEGRALS.C
//  ------------------
//  Version History: 1.0 - 31/12/2015
//
///////////////////////////////////////////////////

#include "stdio.h"
#include "math.h"
#include "scf_globals.h"

void format_integrals(void)
{
	UINT32 i, j, k, l;


	printf("\n\n^^^^Starting the formatting of all integrals into readable matrix form^^^^i\n\n");
	
	//Form the core Hamiltonian - eqn 3.153 pg 141
	h_mat[0][0] = T11 + V11_nucA + V11_nucB;
	h_mat[0][1] = T12 + V12_nucA + V12_nucB;
	h_mat[1][0] = h_mat[0][1];
	h_mat[1][1] = T22 + V22_nucA + V22_nucB;

	//Now for the overlap matrix
	s_mat[0][0] = 1.0;
	s_mat[0][1] = S12;
	s_mat[1][0] = s_mat[0][1];
	s_mat[1][1] = 1.0;

	//We want to use canonical orthogonalisation - eqn 3.238 pg 163 and pg 144
  x_mat[0][0] = 1.0/(sqrt(2.0*(1.0+S12)));
  x_mat[0][1] = 1.0/(sqrt(2.0*(1.0-S12)));
  x_mat[1][0] = x_mat[0][0];
  x_mat[1][1] = -x_mat[0][1];

	//Now transpose this X matrix
	xt_mat[0][0] = x_mat[0][0];
	xt_mat[0][1] = x_mat[1][0];
	xt_mat[1][0] = x_mat[0][1];
	xt_mat[1][1] = x_mat[1][1];
	
	//Matrix of two electron integrals
  tt_mat[0][0][0][0] = V1111;
  tt_mat[1][0][0][0] = V2111;
  tt_mat[0][1][0][0] = V2111;
  tt_mat[0][0][1][0] = V2111;
  tt_mat[0][0][0][1] = V2111;
  tt_mat[1][0][1][0] = V2121;
  tt_mat[0][1][1][0] = V2121;
  tt_mat[1][0][0][1] = V2121;
  tt_mat[0][1][0][1] = V2121;
  tt_mat[1][1][0][0] = V2211;
  tt_mat[0][0][1][1] = V2211;
  tt_mat[1][1][1][0] = V2221;
  tt_mat[1][1][0][1] = V2221;
  tt_mat[1][0][1][1] = V2221;
  tt_mat[0][1][1][1] = V2221;
	tt_mat[1][1][1][1] = V2222;

	if(print_level == PRINT_ALL)
	{
	  //Print the S overlap matrix.
 	 	matout(&s_mat[0][0], 2, 2, sizeof(s_mat), S_MAT_ID);
  	//Print the X transformation matrix.
  	matout(&x_mat[0][0], 2, 2, sizeof(x_mat), X_MAT_ID);
  	//Print the H matrix.
  	matout(&h_mat[0][0], 2, 2, sizeof(h_mat), H_MAT_ID);
	}	

	for(i=0; i<MAX_ELEC; i++)
	{
		for(j=0; j<MAX_ELEC; j++)
		{
			for(k=0; k<MAX_ELEC; k++)
			{
				for(l=0; l<MAX_ELEC; l++)
				{
					printf("(%d %d %d %d) = %.6f\n", i, j, k, l, tt_mat[i][j][k][l]);
				}
			}
		}	
	}

	printf("\n^^^^Completed the formatting of integrals into readable format^^^^\n\n");
}


