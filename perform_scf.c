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
//  PERFORM_SCF.C
//  -------------
//  Version History: 1.0 - 31/12/2015
//
///////////////////////////////////////////////////

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "scf_globals.h"

//Convergence criteria and maximum number of SCF iterations
double conv_crit = 0.0001;
UINT32 maxiter = 25;
BOOLEAN converged = FALSE;

void formg(void);
void mult_mat (double *matA, double *matB, UINT32 m,  double *result_mat);
void diag(void);

//This converts the MATRIX_ID to a string.
//Always make this lookup table match MATRIX_ID
char mat_id_2_str[LAST_MAT_ID+1][MAX_STR] =
{
  "H",
	"S",
  "G",
  "P",
  "E",
  "C",
  "F",
  "Old P",
  "F Prime",
  "C Prime",
	"X",
	"XT",
	"MULLIKEN",
	"TEMP_2D",

  //ALWAYS MAKE THIS THE LAST LINE
  "END OF MATRIX LIST"
};


void perform_scf(void)
{
  UINT32 iter;
  UINT32 i, j, k;
	double conv = 100.0; //Pick a large number to ensure above convergence threshold
	double delta;


	printf("\n\n----Performing SCF----\n\n");
	
	//Use the core Hamilton as the starting guess for the density matrix
	//This is SCF step 4 page 146
  for (i=0; i<MAX_ELEC; i++)
  {
		for (j=0; j<MAX_ELEC; j++)
		{
			//Initial guess of core Hamiltonian means the density matrix is all zeros
			dens_mat[i][j] = 0.0;			
		}
	}
	
	iter = 1;

  while((iter <= maxiter) && (!converged))
	{

		if(print_level == PRINT_ALL)
		{
      printf("Starting Iteration %d.....\n\n", iter);

		  //Print the density matrix P.
	 	  matout(&dens_mat[0][0], 2, 2, sizeof(dens_mat), P_MAT_ID);
		}
	
		//For the 2 electron G matrix as per step 5 page 146
		formg();			

	  if(print_level == PRINT_ALL)
 	 	{
 	   	//Print the density matrix P.
 	   	matout(&g_mat[0][0], 2, 2, sizeof(g_mat), G_MAT_ID);
      matout(&h_mat[0][0], 2, 2, sizeof(h_mat), H_MAT_ID);
		}


		//Add the G matrix to the core Hamiltonian to get the Fock matrix
		//Eqn 3.154 page 141 and step 6 page 146
		for (i=0; i<MAX_ELEC; i++)
		{
			for (j=0; j<MAX_ELEC; j++)
			{
				f_mat[i][j] = h_mat[i][j] + g_mat[i][j];
			}
		}

		//Calculate the energy now.
		//Eqn 3.184 pg 150
		elec_energy = 0.0;
		for(i=0; i<MAX_ELEC; i++)
		{
			for (j=0; j<MAX_ELEC; j++)
			{
				//Eqn 3.184 pg 150
				elec_energy = elec_energy + 0.5*dens_mat[i][j]*(h_mat[i][j] + f_mat[i][j]);
			}
		}
				
    if(print_level == PRINT_ALL)
    {
      //Print the density matrix P.
      matout(&f_mat[0][0], 2, 2, sizeof(f_mat), F_MAT_ID);

			printf("\n------------------------------\n");
			printf("Intermediate Electronic Energy = %.6f\n", elec_energy);
      printf("-------------------------------\n\n");
    }		

		//Transform the Fock matrix. FPRIME = XT*F*X (step 7 pg 146)
		mult_mat(&f_mat[0][0], &x_mat[0][0], 2, &temp_2d_mat[0][0]);		
    mult_mat(&xt_mat[0][0], &temp_2d_mat[0][0], 2, &fprime_mat[0][0]);

		//Diagonalise the Transformed Fock Matrix
		//Step 8 pg 146
		diag();
	
		//Now back-transform the eigenvector matrix to get the matrix C
		//Step 9 pg 146
    mult_mat(&x_mat[0][0], &cprime_mat[0][0], 2, &c_mat[0][0]);

		//Form new density matrix
		//Step 10 pg 146
		for (i=0; i<MAX_ELEC; i++)
		{
			for(j=0; j<MAX_ELEC; j++)
			{
				//Save current density matrix before creating a new one.
				olddens_mat[i][j] = dens_mat[i][j];
				dens_mat[i][j] = 0.0;

				//Not sure why the following loop stops at k=1. We only have one occupied MO
				//though as this is a restricted calculation with 2 electrons so that might be it.
				//eqn 3.145 pg 139
				for (k=0; k<1; k++)
				{
					dens_mat[i][j] = dens_mat[i][j] + 2.0*c_mat[i][k]*c_mat[j][k];
				}
			}
		}
		
    if(print_level == PRINT_ALL)
    {
			//Print the Fprime matrix
      matout(&fprime_mat[0][0], 2, 2, sizeof(fprime_mat), FPRIME_MAT_ID);
      matout(&cprime_mat[0][0], 2, 2, sizeof(cprime_mat), CPRIME_MAT_ID);
      matout(&energy_mat[0][0], 2, 2, sizeof(energy_mat), E_MAT_ID);
      matout(&c_mat[0][0], 2, 2, sizeof(c_mat), C_MAT_ID);
      matout(&dens_mat[0][0], 2, 2, sizeof(dens_mat), P_MAT_ID);
    }

		//Now for convergence checking get the delta between the new density matrix and the previous one.
		//Use least squares method - step 11 pg 146.
		delta = 0.0;
		for (i=0;i<MAX_ELEC; i++)
		{
			for(j=0; j<MAX_ELEC; j++)
			{
				delta = delta + pow(dens_mat[i][j] - olddens_mat[i][j], 2);
			}
		}
		delta = sqrt(delta/4);
		printf("DeltaP measured as %.6f\n\n", delta);
		
		//Now check for convergence
		if (delta < conv_crit)
		{
			//We have reached convergence. Set flag to get out of here. 
			converged = TRUE;
		}
		else
		{	
			//Not converged yet. Further iterations needed.
			iter++;
		}
	}

	//Did we get here because we converged or because we ran out of iterations?
	if(converged == TRUE)
	{
		//It's good, we converged.
    printf("Convergence achieved in %d steps\n", iter);
	}
	else if (iter > maxiter)
	{
		//We ran out of iterations
		printf("Number of iterations %d exceeds maximum of %d\n", iter, maxiter);
	}
	else
	{
		printf("Fucked something up because we got here without exceeding iterations or converging!!!\n");
	}

	//Add the nuclear energy contribution to the electronic energy we have calculated.
	//eqn 2.14 pg 44
	total_energy = elec_energy + (za*zb)/bond_len;

	//Calculate Mulliken Populations
	//eqn 3.305 pg 203
  mult_mat(&dens_mat[0][0], &s_mat[0][0], 2, &mulliken_mat[0][0]);
	
	//Now print out the results and the final matrices
	printf("\n++++++++++FINAL RESULTS++++++++++\n\n");
  matout(&g_mat[0][0], 2, 2, sizeof(g_mat), G_MAT_ID);
  matout(&f_mat[0][0], 2, 2, sizeof(f_mat), F_MAT_ID);
	matout(&energy_mat[0][0], 2, 2, sizeof(energy_mat), E_MAT_ID);
  matout(&c_mat[0][0], 2, 2, sizeof(c_mat), C_MAT_ID);
  matout(&dens_mat[0][0], 2, 2, sizeof(dens_mat), P_MAT_ID);
  printf("\n\n***************************************\n\n");
  printf("Final Electronic Energy = %.6f\n", elec_energy);
  printf("Final Total Energy = %.6f\n\n", total_energy);
  matout(&mulliken_mat[0][0], 2, 2, sizeof(mulliken_mat), MULLIKEN_MAT_ID);
	printf("\n***************************************\n\n");
  printf("\n+++++++++++++++++++++++++++++++++\n");

	printf("\n----End of SCF Procedure----\n");
}


void formg(void)
{
	UINT32 i, j, k, l;

	for(i=0; i<MAX_ELEC; i++)
	{
		for(j=0; j<MAX_ELEC; j++)
		{
			g_mat[i][j] = 0.0;

			for(k=0; k<MAX_ELEC; k++)
			{
				for(l=0; l<MAX_ELEC; l++)
				{
					//Eqn 3.154 page 141
					g_mat[i][j] = g_mat[i][j] + (dens_mat[k][l]*(tt_mat[i][j][k][l] - (0.5*tt_mat[i][l][k][j])));
				}
			}
		}
	}
}


void mult_mat (double *matA, double *matB, UINT32 m,  double *result_mat)
{
	//Multiplies two m x m square matrices together.
	UINT32 i, j, k;

	for (i=0; i<m; i++)
	{
		for (j=0; j<m; j++)
		{
			*(result_mat + i*m + j) = 0.000000;			

			for (k=0; k<m; k++)
			{
				*(result_mat + i*m + j) = *(result_mat + i*m + j) + (*(matA + i*m + k) * (*(matB + k*m + j)));
			}
		}
	}
}

void diag(void)
{
	double theta;
	double temp;	

	
	//Diagonalises the Transformed Fock matrix.
	//In this program we are using the 2 x 2 matrix solution
	//from chapter 1.
	if(fabs(fprime_mat[0][0] - fprime_mat[1][1]) < pow(10, -20))
	{
		//Symmetry determined solution. We have to do this because eqn 1.105 pg 20
		//divides by the difference of these two matrix elements. Must avoid divide
		//by zero. For zero difference between them we know that is simply an asymptote
		//at 90 degrees or PI/4 so simply set that angle directly to stop the program
		//crashing.
		theta = PI/4;
	}
	else
	{
		//eqn 1.105 pg 20		
		theta = 0.5*atan(2.0*fprime_mat[0][1]/(fprime_mat[0][0]-fprime_mat[1][1]));
	}

	//Fill in the eigenvector matrix
	//eqn 1.104 pg 20
	cprime_mat[0][0] = cos(theta);  
	cprime_mat[0][1] = sin(theta);
  cprime_mat[1][0] = sin(theta);
  cprime_mat[1][1] = -cos(theta);

	//Finally we can now calculate the energy matrix directly.
	//eqn 1.106a and 1.106b pg 21. Off diags are zero.
  energy_mat[0][0] = fprime_mat[0][0]*cos(theta)*cos(theta) + fprime_mat[1][1]*sin(theta)*sin(theta) 
												+ fprime_mat[0][1]*sin(2*theta);
  energy_mat[0][1] = 0.0;
	energy_mat[1][0] = 0.0;
  energy_mat[1][1] = fprime_mat[1][1]*cos(theta)*cos(theta) + fprime_mat[0][0]*sin(theta)*sin(theta) 
												- fprime_mat[0][1]*sin(2*theta);

	//Now sort the eigenvalues and eigenvectors in ascending order if they are out of order.
	if (energy_mat[1][1] < energy_mat[0][0])
	{
		//In wrong order. Swap eigenvalues and eigenvectors
		temp = energy_mat[1][1];
		energy_mat[1][1] = energy_mat[0][0];
		energy_mat[0][0] = temp;

		temp = cprime_mat[0][1];
		cprime_mat[0][1] = cprime_mat[0][0];
		cprime_mat[0][0] = temp;
		temp = cprime_mat[1][1];
		cprime_mat[1][1] = cprime_mat[1][0];
		cprime_mat[1][0] = temp;
	}
}



void matout(double *matrix, UINT32 m, UINT32 n, UINT32 mat_size, MATRIX_ID mat_id)
{
	//Print out an m x n matrix

	UINT32 i, j;

	
	//First ensure we have passed the right arguments or we risk overrunning matrix.
  if(mat_size != (m*n*sizeof(double)))
	{
		printf("Error reading matrix. Size is %d, m is %d, n is %d, mat_id is %s \n", mat_size, m, n, mat_id_2_str[mat_id]);
	}

	printf("The %s Matrix\n", mat_id_2_str[mat_id]);
  for(j=0; j<n; j++)
	{
		printf("           %d", j+1);
	}
	printf("\n");

	for(i=0; i<m; i++)
	{
		printf("%d", i);		

		for(j=0; j<n; j++)
		{
			printf("%14.6f ", *(matrix + i*n + j));
		}

		printf("\n\n");
	}
}
