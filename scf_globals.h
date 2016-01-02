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
//  SCF_GLOBALS.H
//  -------------
//  Version History: 1.0 - 31/12/2015
//
///////////////////////////////////////////////////

typedef unsigned int UINT32;

typedef enum
{
	FALSE = 0,
	TRUE
} BOOLEAN;

//Print levels
typedef enum
{
	PRINT_NONE = 0,						//No printing
	PRINT_CONV_RESULTS_ONLY,	//Print only the final converged results
	PRINT_ALL									//Print everything.
} PRINT_LEVELS;

PRINT_LEVELS print_level;

//Select the n in STO-nG basis set
UINT32 sto_ng; 
 
//Matrix Identifiers
//We use this for printing so we know what text to add to the
//matrix printout
typedef enum
{
  H_MAT_ID = 0,
  S_MAT_ID,
	G_MAT_ID,
	P_MAT_ID,
	E_MAT_ID,
	C_MAT_ID,
	F_MAT_ID,
	OLDP_MAT_ID,
	FPRIME_MAT_ID,
	CPRIME_MAT_ID,
	X_MAT_ID,
	XT_MAT_ID,
	MULLIKEN_MAT_ID,
	TEMP_2D_MAT_ID,

	//ALWAYS MAKE THIS THE LAST LINE
	LAST_MAT_ID
  
} MATRIX_ID;

#define MAX_STR 132

char mat_id_2_str[LAST_MAT_ID+1][MAX_STR]; 

//Only supporting up to STO-3G at the moment
#define MAX_STO 3

//Define PI to 6 dp.
#define PI 3.141593

//Define the maximum supported number of electrons
#define MAX_ELEC 2

double bond_len; //bond length
double bond_len_sq; //bond length squared
double zetaA; //slater orbital exponent for atom A
double zetaB; //slater orbital exponent for atom B
double za; //nuclear charge of atom A
double zb; //nuclear charge of atom B
double elec_energy; //calculated electronic energy
double total_energy; //calculated electronic energy + nuclear energy contribution

//Store the integrals
double S12, T11, T12, T22, V11_nucA, V12_nucA, V22_nucA, V11_nucB, V12_nucB, V22_nucB;
double V1111, V2111, V2121, V2211, V2221, V2222;

//Global matrices
double s_mat[MAX_ELEC][MAX_ELEC], x_mat[MAX_ELEC][MAX_ELEC], xt_mat[MAX_ELEC][MAX_ELEC]; 
double h_mat[MAX_ELEC][MAX_ELEC], f_mat[MAX_ELEC][MAX_ELEC], g_mat[MAX_ELEC][MAX_ELEC];
double c_mat[MAX_ELEC][MAX_ELEC], fprime_mat[MAX_ELEC][MAX_ELEC], cprime_mat[MAX_ELEC][MAX_ELEC];
double dens_mat[MAX_ELEC][MAX_ELEC], olddens_mat[MAX_ELEC][MAX_ELEC], energy_mat[MAX_ELEC][MAX_ELEC];
double mulliken_mat[MAX_ELEC][MAX_ELEC], temp_2d_mat[MAX_ELEC][MAX_ELEC];
double tt_mat[MAX_ELEC][MAX_ELEC][MAX_ELEC][MAX_ELEC];

//function declarations
void calc_integrals(void);
void format_integrals(void);
void perform_scf(void);
void matout(double *matrix, UINT32 m, UINT32 n, UINT32 mat_size, MATRIX_ID mat_id);

