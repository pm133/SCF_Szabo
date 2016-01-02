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
//  CALC_INTEGRALS.C
//  ----------------
//  Version History: 1.0 - 31/12/2015
//
///////////////////////////////////////////////////

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "scf_globals.h"

//This function is where the exponents of the Gaussians primitives are scaled
//and the contraction coefficients are setup to approximate a 1s Slater orbital.
void prepare_basis_fns(void);

void initialise_integrals(void);

//Entry point for evaluating one electron integrals
void eval_one_elec_integrals(void);
double overlap(double nucA_expon, double nucB_expon, double r2);
double kinetic_energy(double nucA_expon, double nucB_expon, double r2);
double one_elec_nuc_attr(double nucA_expon, double nucB_expon, double rab2, double rcp2, double zc);

//Calculates erf function
double FO(double arg);

//Entry point for evaluating two electron integrals
void eval_two_elec_integrals(void);
//Calculates the 2 electron repulsions
double two_elec_repulsion(double nucA_expon, double nucB_expon, double nucC_expon, double nucD_expon, double rab2, double rcd2, double rpq2);

//Normalised coefficients for unscaled exponents (zeta = 0) of STO-NG basis. 
//Each row is for a specific value of N.
//STO-1G, STO-2G and STO-3G supported. Note that STO-1G and STO-2G have zeroes
//where no coefficients exist.
//These coeffs are found on page 157 eqns 3.219 to 3.221.
double coef_sto_ng_basis [MAX_STO][MAX_STO] =
{
	1.000000, 0.000000, 0.000000,
	0.678914, 0.430129, 0.000000,
	0.444635, 0.535328, 0.154329
};

//Unscaled (zeta = 1.0) exponents of STO-NG basis. Each row is for a specific value of N.
//STO-1G, STO-2G and STO-3G supported. Note that STO-1G and STO-2G have zeroes
//where no coefficients exist.
//These coeffs are found on page 157 eqns 3.219 to 3.221.
double exp_sto_ng_basis [MAX_STO][MAX_STO] =
{
	0.270950, 0.000000, 0.000000,
	0.151623, 0.851819, 0.000000,
	0.109818, 0.405771, 2.227660
};

//This set of arrays store the scaled coefficients for the nuclei depending on their zeta value.
double scaled_exp_nucA[MAX_STO];
double scaled_exp_nucB[MAX_STO];

//Now these arrays store the re-normalised coefficients based on the newly scaled exponents.
double renorm_coefs_nucA[MAX_STO];
double renorm_coefs_nucB[MAX_STO];

//Define some distances and their squared values for integrals
//R_ap is the distance between centre A and centre P etc.
double R_ap, R_ap_sq, R_bp, R_bp_sq, R_aq, R_aq_sq, R_bq, R_bq_sq, R_pq, R_pq_sq;  

//This is the function responsible for performing the integrals.
void calc_integrals(void)
{
	printf("\n\n++++Starting the calcuation of all integrals++++\n");

	//The first job is to sort out the exponents and the contraction coefficients of the basis functions.
	//First scale the exponents for the corresponding tau value for each nuclear centre.
	//Then need to re-normalise the coefficients based on these exponents.
	prepare_basis_fns();

  bond_len_sq = bond_len*bond_len;

	//Default integral elements to zero
	initialise_integrals();

	//Now evaluate the various integrals
	eval_one_elec_integrals();
	eval_two_elec_integrals();	

	if(print_level == PRINT_ALL)
	{
		printf("\n\n");
    printf("Bond Length = %.6f\n", bond_len);
		printf("Bond Length Squared = %.6f\n", bond_len_sq);
    printf("Zeta_A = %.6f\n", zetaA);
    printf("Zeta_B = %.6f\n", zetaB);
    printf("S12 = %.6f\n", S12);
    printf("T11 = %.6f\n", T11);
    printf("T12 = %.6f\n", T12);
    printf("T22 = %.6f\n", T22);
    printf("V11_nucA = %.6f\n", V11_nucA);
    printf("V12_nucA = %.6f\n", V12_nucA);
    printf("V22_nucA = %.6f\n", V22_nucA);
    printf("V11_nucB = %.6f\n", V11_nucB);
    printf("V12_nucB = %.6f\n", V12_nucB);
    printf("V22_nucB = %.6f\n", V22_nucB);
    printf("V1111 = %.6f\n", V1111);
    printf("V2111 = %.6f\n", V2111);
    printf("V2121 = %.6f\n", V2121);
    printf("V2211 = %.6f\n", V2211);
    printf("V2221 = %.6f\n", V2221);
    printf("V2222 = %.6f\n", V2222);
	}

	printf("\n++++Completed calculation of all integrals++++\n");
}


void prepare_basis_fns(void)
{
	//Create the contraction coefficients and exponents to approximate a normalised 1s Slater orbital
	//with exponent 1.0 using normalised 1s Gaussian type primitives
	
	UINT32 i;

	if(print_level = PRINT_ALL)
	{
		for(i=0; i<sto_ng; i++)
		{
			printf("coef[%d][%d] = %.6f; ", sto_ng-1, i, coef_sto_ng_basis[sto_ng-1][i]);		
			printf("        ");
		}
		
		printf("\n\n");

		for(i=0; i<sto_ng; i++)
		{
      printf("exp[%d][%d] = %.6f; ", sto_ng-1, i, exp_sto_ng_basis[sto_ng-1][i]);
      printf("        ");
		}	
   	printf("\n\n");
	}

	//Now scale exponents for the zeta values appropriate for each nuclear centre.
	//The equations for this come from page 158 eqn 3.224 for the scaling and 
	//page 153 eqn 3.203 for the renormalisation.
	for(i=0; i<sto_ng; i++)
	{
		scaled_exp_nucA[i] = (zetaA*zetaA)*exp_sto_ng_basis[sto_ng-1][i];
   	scaled_exp_nucB[i] = (zetaB*zetaB)*exp_sto_ng_basis[sto_ng-1][i];
    renorm_coefs_nucA[i] = pow((2.0*scaled_exp_nucA[i]/PI), 0.75)*coef_sto_ng_basis[sto_ng-1][i]; 
    renorm_coefs_nucB[i] = pow((2.0*scaled_exp_nucB[i]/PI), 0.75)*coef_sto_ng_basis[sto_ng-1][i];
	}


  if(print_level = PRINT_ALL)
  {
    for(i=0; i<sto_ng; i++)
    {
      printf("scaled_exp_nucA[%d] = %.6f; ", i, scaled_exp_nucA[i]);
  	}

    printf("\n\n");
    for(i=0; i<sto_ng; i++)
    {
      printf("scaled_exp_nucB[%d] = %.6f; ", i, scaled_exp_nucB[i]);
    }

    printf("\n\n");
    for(i=0; i<sto_ng; i++)
    {
      printf("renorm_coefs_nucA[%d] = %.6f; ", i, renorm_coefs_nucA[i]);
    }

    printf("\n\n");
    for(i=0; i<sto_ng; i++)
    {
      printf("renorm_coefs_nucB[%d] = %.6f; ", i, renorm_coefs_nucB[i]);
    }

    printf("\n\n");
  }

}	


void initialise_integrals(void)
{
	//Only 2 electrons in this system so easier and more efficient to hold
	//a separate variable for each array element.

	//There is only one overlap element as there are only two electrons
	S12 = 0.0;

	//Only storing triangle as matrix real and symmetric
	T11 = 0.0;
  T12 = 0.0;
  T22 = 0.0;

  //Only storing triangle as matrix real and symmetric
  V11_nucA = 0.0;
  V12_nucA = 0.0;
  V22_nucA = 0.0;

  //Only storing triangle as matrix real and symmetric
  V11_nucB = 0.0;
  V12_nucB = 0.0;
  V22_nucB = 0.0;

  //Only storing triangle as matrix real and symmetric
  V1111 = 0.0;
  V2111 = 0.0;
  V2121 = 0.0;
  V2211 = 0.0;
  V2221 = 0.0;
  V2222 = 0.0;
}

void eval_one_elec_integrals(void)
{
	UINT32 i, j;

	//Perform calculations of the one electron integrals here.
	//Nucleus A is the first atom, nucleus B is the second atom.
	//Origin on nucleus A.
	//V12_nucA is the off diagonal nuclear attraction to nucleus A etc.
	for(i=0; i<sto_ng; i++)
	{
		for(j=0; j<sto_ng; j++)
		{
			R_ap = (scaled_exp_nucB[j]*bond_len)/(scaled_exp_nucA[i]+scaled_exp_nucB[j]);
			R_ap_sq = R_ap*R_ap;
			R_bp_sq = pow((bond_len - R_ap), 2.0);

			//This is equation A.9 from page 412
			S12 = S12 + renorm_coefs_nucA[i] * renorm_coefs_nucB[j] * 
						overlap(scaled_exp_nucA[i], scaled_exp_nucB[j], bond_len_sq);

 			//This is equation A.11 from page 412
			T11 = T11 + renorm_coefs_nucA[i] * renorm_coefs_nucA[j] * 
						kinetic_energy(scaled_exp_nucA[i], scaled_exp_nucA[j], 0.0);

      //This is equation A.11 from page 412
      T12 = T12 + renorm_coefs_nucA[i] * renorm_coefs_nucB[j] *
					 	kinetic_energy(scaled_exp_nucA[i], scaled_exp_nucB[j], bond_len_sq);

      //This is equation A.11 from page 412
      T22 = T22 + renorm_coefs_nucB[i] * renorm_coefs_nucB[j] *
 						kinetic_energy(scaled_exp_nucB[i], scaled_exp_nucB[j], 0.0);

			//Now we need the attractions between nuclei and electrons
			//This is equation A.32 and A.33 on page 415
			V11_nucA = V11_nucA + renorm_coefs_nucA[i] * renorm_coefs_nucA[j] * 
						one_elec_nuc_attr(scaled_exp_nucA[i], scaled_exp_nucA[j], 0.0, 0.0, za);

      //This is equation A.32 and A.33 on page 415
      V12_nucA = V12_nucA + renorm_coefs_nucA[i] * renorm_coefs_nucB[j] * 
            one_elec_nuc_attr(scaled_exp_nucA[i], scaled_exp_nucB[j], bond_len_sq, R_ap_sq, za);

      //This is equation A.32 and A.33 on page 415
      V22_nucA = V22_nucA + renorm_coefs_nucB[i] * renorm_coefs_nucB[j] * 
            one_elec_nuc_attr(scaled_exp_nucB[i], scaled_exp_nucB[j], 0.0, bond_len_sq, za);

      //This is equation A.32 and A.33 on page 415
      V11_nucB = V11_nucB + renorm_coefs_nucA[i] * renorm_coefs_nucA[j] * 
            one_elec_nuc_attr(scaled_exp_nucA[i], scaled_exp_nucA[j], 0.0, bond_len_sq, zb);

      //This is equation A.32 and A.33 on page 415
      V12_nucB = V12_nucB + renorm_coefs_nucA[i] * renorm_coefs_nucB[j] * 
            one_elec_nuc_attr(scaled_exp_nucA[i], scaled_exp_nucB[j], bond_len_sq, R_bp_sq, zb);
        
      //This is equation A.32 and A.33 on page 415
      V22_nucB = V22_nucB + renorm_coefs_nucB[i] * renorm_coefs_nucB[j] * 
            one_elec_nuc_attr(scaled_exp_nucB[i], scaled_exp_nucB[j], 0.0, 0.0, zb);
		}
	}
}

double overlap(double nucA_expon, double nucB_expon, double r2)
{
	//Part of eqn A.9 on page 412
	return pow((PI/(nucA_expon + nucB_expon)), 1.5) * 
				exp(-nucA_expon * nucB_expon * 
					r2/(nucA_expon + nucB_expon));	
}

double kinetic_energy(double nucA_expon, double nucB_expon, double r2)
{
	//Part of eqn A.11 on page 411
	return (nucA_expon * nucB_expon/(nucA_expon + nucB_expon)) *
          	 (3.0 - (2.0 * nucA_expon * nucB_expon * r2/(nucA_expon + nucB_expon))) *
            	 pow(PI/(nucA_expon + nucB_expon), 1.5) *
                 	exp(-nucA_expon * nucB_expon * r2/(nucA_expon + nucB_expon));

}

double one_elec_nuc_attr(double nucA_expon, double nucB_expon, double rab2, double rcp2, double zc)
{
 	//This is equation A.32 and A.33 on page 415
 	return 2.0 * (PI/(nucA_expon + nucB_expon)) * 
				FO((nucA_expon + nucB_expon) * rcp2) *
					exp((-1)*nucA_expon * nucB_expon * rab2/(nucA_expon + nucB_expon)) * 
						((-1)*zc);
}

double FO(double arg)
{
	//Here 'arg' is just acomplicated expression that we're trying to the erf function.
	//Calling erf with values under 1*10^(-6)
	if (arg >= 0.000001)
	{
//		printf("FO-A\n");
		return 0.5 * sqrt(PI/arg) * erf(sqrt(arg));
	}
	else
	{
//        printf("FO-B\n");
		return (1.0 - (arg/(double)3));
	}
}

void eval_two_elec_integrals(void)
{
	UINT32 i, j, k, l;

	//Now evaluate the two electron integrals
	//This quadruple loop is at the heart of the problem!!!
	for(i=0; i<sto_ng; i++)
	{
		for(j=0; j<sto_ng; j++)
		{
			for(k=0; k<sto_ng; k++)
			{
				for(l=0; l<sto_ng; l++)
				{
					//Get the various bond lengths and square them
					R_ap = scaled_exp_nucB[i] * bond_len / (scaled_exp_nucB[i] + scaled_exp_nucA[j]);
					R_bp = bond_len - R_ap;
					R_aq = scaled_exp_nucB[k] * bond_len / (scaled_exp_nucB[k] + scaled_exp_nucA[l]);
					R_bq = bond_len - R_aq;
					R_pq = R_ap - R_aq;
					R_ap_sq = R_ap * R_ap;
					R_bp_sq = R_bp * R_bp;
					R_aq_sq = R_aq * R_aq;
					R_bq_sq = R_bq * R_bq;
					R_pq_sq = R_pq * R_pq;

					//These formulae are equation A.32 on page 415 and equation A.41 page 416
					V1111 = V1111 + renorm_coefs_nucA[i] * renorm_coefs_nucA[j] * renorm_coefs_nucA[k] * renorm_coefs_nucA[l] *
																two_elec_repulsion(scaled_exp_nucA[i], scaled_exp_nucA[j], scaled_exp_nucA[k], scaled_exp_nucA[l], 
																	0.0, 0.0, 0.0);

					V2111 = V2111 + renorm_coefs_nucB[i] * renorm_coefs_nucA[j] * renorm_coefs_nucA[k] * renorm_coefs_nucA[l] *
                                two_elec_repulsion(scaled_exp_nucB[i], scaled_exp_nucA[j], scaled_exp_nucA[k], scaled_exp_nucA[l], 
                                  bond_len_sq, 0.0, R_ap_sq);

					V2121 = V2121 + renorm_coefs_nucB[i] * renorm_coefs_nucA[j] * renorm_coefs_nucB[k] * renorm_coefs_nucA[l] *
                                two_elec_repulsion(scaled_exp_nucB[i], scaled_exp_nucA[j], scaled_exp_nucB[k], scaled_exp_nucA[l], 
                                  bond_len_sq, bond_len_sq, R_pq_sq);

					V2211 = V2211 + renorm_coefs_nucB[i] * renorm_coefs_nucB[j] * renorm_coefs_nucA[k] * renorm_coefs_nucA[l] *
                                two_elec_repulsion(scaled_exp_nucB[i], scaled_exp_nucB[j], scaled_exp_nucA[k], scaled_exp_nucA[l],
																	0.0, 0.0, bond_len_sq); 

					V2221 = V2221 + renorm_coefs_nucB[i] * renorm_coefs_nucB[j] * renorm_coefs_nucB[k] * renorm_coefs_nucA[l] *
                                two_elec_repulsion(scaled_exp_nucB[i], scaled_exp_nucB[j], scaled_exp_nucB[k], scaled_exp_nucA[l], 
                                  0.0, bond_len_sq, R_bq_sq);

					V2222 = V2222 + renorm_coefs_nucB[i] * renorm_coefs_nucB[j] * renorm_coefs_nucB[k] * renorm_coefs_nucB[l] *
                                two_elec_repulsion(scaled_exp_nucB[i], scaled_exp_nucB[j], scaled_exp_nucB[k], scaled_exp_nucB[l], 
                                  0.0, 0.0, 0.0);
				}
			}
		}
	}
}


double two_elec_repulsion(double nucA_expon, double nucB_expon, double nucC_expon, double nucD_expon, double rab2, double rcd2, double rpq2)
{
	return 2.0 * pow(PI,2.5)/((nucA_expon + nucB_expon) *
				(nucC_expon + nucD_expon) * sqrt (nucA_expon + nucB_expon + nucC_expon + nucD_expon)) * 
					FO((nucA_expon + nucB_expon) * (nucC_expon + nucD_expon) * rpq2/(nucA_expon + nucB_expon + nucC_expon + nucD_expon)) *
						exp((-nucA_expon*nucB_expon*rab2/(nucA_expon + nucB_expon)) - (nucC_expon*nucD_expon*rcd2/(nucC_expon+nucD_expon)));	
}


