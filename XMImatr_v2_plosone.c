/********************************************************************************/

/*
 * If you use this toolbox, *please cite* the following references:
 *
 * C. Rummel, M. Müller, G. Baier, F. Amor, K. Schindler
 * Analyzing spatio-temporal patterns of genuine cross-correlations
 * J. Neurosci. Meth. 191, 94-100 (2010)
 *
 * C. Rummel, E. Abela, M. Müller, M. Hauf, O. Scheidegger, R. Wiest, K. Schindler
 * Uniform approach to linear and nonlinear interrelation patterns in multivariate time series
 * Phys. Rev. E83, 066215 (2011)
 *
 * C. Rummel, E. Abela, R.G. Andrzejak, M. Hauf, C. Pollo, M. Müller, C. Weisstanner, R. Wiest, K. Schindler
 * Resected brain tissue, seizure onset zone and quantitative EEG measures:
 * Towards prediction of post-surgical seizure control
 * PLoS ONE accepted (2015)
 */

/*
 * This is the main routine used for reading of data, calculation of surrogates
 * and matrices, statistics and multiple comparison correction and
 * writing of output.
 */

/*
 * In addition to the distributed source code, the following software is required:
 * - the MILCA package (can be downloaded at
 *   https://www.ucl.ac.uk/ion/departments/sobell/Research/RLemon/MILCA/MILCA,
 *   status 2015/10/12) for calculation of mutual information
 * - Numerical Recipes in C (license required).
 */

/*
 * several INPUT  files are required in the directory ../inputs (must be provided)
 * several OUTPUT files are written to  the directory ../data   (must be provided)
 */

/*
 * written by Christian Rummel, 2010/08/17, crummel@web.de
*/

/********************************************************************************/


// include files

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <time.h>
#include <sys/types.h>



int main ()
{
// data types

FILE *fout_TEST;

	FILE *fin_calc;						// input file for calculation parameters
	FILE *fin_series;					// input file for time series parameters and data

	FILE *fout_matrix_TCS;				// output file for median of matrix elements
	FILE *fout_matrix_CCS;				// output file for median of matrix elements
	FILE *fout_matrix_RCS;				// output file for median of matrix elements
	FILE *fout_matrix_multCCS;			// output file for median of matrix elements
	FILE *fout_matrix_multRCS;			// output file for median of matrix elements

	FILE *fout_matrix_TMI;				// output file for median of matrix elements
	FILE *fout_matrix_CMI;				// output file for median of matrix elements
	FILE *fout_matrix_RMI;				// output file for median of matrix elements
	FILE *fout_matrix_multCMI;			// output file for median of matrix elements
	FILE *fout_matrix_multRMI;			// output file for median of matrix elements

	char name_fout_matrix_TCS    [90] = "../data/";
	char name_fout_matrix_CCS    [90] = "../data/";
	char name_fout_matrix_RCS    [90] = "../data/";
	char name_fout_matrix_multCCS[90] = "../data/";
	char name_fout_matrix_multRCS[90] = "../data/";

	char name_fout_matrix_TMI    [90] = "../data/";
	char name_fout_matrix_CMI    [90] = "../data/";
	char name_fout_matrix_RMI    [90] = "../data/";
	char name_fout_matrix_multCMI[90] = "../data/";
	char name_fout_matrix_multRMI[90] = "../data/";

	int i, ii, j, k, n, nn;				// indices

	int neighbors;						// number of nearest neighbors used for mutual information estimation
										// if ==0 only cross-correlation based interrelation is estimated

	int T_surr;							// length of surrogate window, must be 2^n
	int T;								// length of analysis window,  must be 2^n
	int M;								// number of time series/channels
	int N_ens;							// ensemble size for statistics
	int N_surr;							// number of surrogates to be generated
	int step_surr, step;				// step width for output in sample points

	unsigned char *line;				// line for successive reading of input files

	int sw_time_in, sw_time_out;		// column for time, no time if ==0
	int sw_slopes;						// switch for analysis of slope time series

	int *sel;							// vector for selection of time series included in analysis

	double f_samp;						// sampling frequency
	int columns;						// number of columns
	double tstart, tstop;				// starting and termination point for analysis

	long seed;							// seed for random numbers
	time_t t1;							// time to set seed for random numbers
	double siglev;						// level of significance
	int sw_corr;						// switch for correction

	double t;							// time argument
	int l;								// discrete time argument

	double readdata;					// dummy for reading of data
	double *Xl, *Xl_old, **X;			// vector and matrix for original time series
	double **X_read, **X_surr;			// matrix for reading and generating surrogates from original time series
	double ave, dev;					// moments for normalization
	double **Xtilde;					// matrix for normalized time series

	int **iter_uni;						// matrix for necessary iterations
	double **accur_uni;					// matrix for reached accuracy
	double ***X_IFS, **X_IFSsing;		// cube and matrix for Independent Fourier Based Surrogate time series
	double **C_IFS, **MI_IFS;			// interrelation matrices of Independent Fourier Based Surrogate time series

	int *iter_mult;						// vector for necessary iterations
	double **accur_mult;				// matrix for reached accuracy
	double ***X_MFS, **X_MFSsing;		// cube and matrix for Multivariate Fourier Based Surrogate time series
	double **C_MFS, **MI_MFS;			// interrelation matrices of Multivariate Fourier Based Surrogate time series

	double **C, **C_surr;				// equal time correlation matrix
	double ***all_C, ***all_C_uni, ***all_C_mult;// cubes for equal time correlation matrices
	double ***all_Cabs, ***all_Cabs_uni, ***all_Cabs_mult;// cubes for equal time correlation matrix
	double *sort_C, *sort_C_surr;		// vectors for sorted equal time correlation matrix elements
	double *sort_Cabs, *sort_Cabs_uni, *sort_Cabs_mult;	// vectors for sorted equal time correlation matrix elements

	double **MI, **MI_surr;				// mutual information
	double ***all_MI, ***all_MI_uni, ***all_MI_mult;// cubes for mutual information matrices
	double ***all_MIabs, ***all_MIabs_uni, ***all_MIabs_mult;// cubes for mutual information matrices
	double *sort_MI, *sort_MI_surr;		// vectors for sorted equal time correlation matrix elements
	double *sort_MIabs, *sort_MIabs_uni, *sort_MIabs_mult;	// vectors for sorted equal time correlation matrix elements

	double *vec_med;					// vector for definition of median
	double **C_med, **C_uni_med, **C_mult_med;// equal time correlation matrix
	double **Cabs_med, **Cabs_uni_med, **Cabs_mult_med;// equal time correlation matrix

	double **MI_med, **MI_uni_med, **MI_mult_med;// equal time correlation matrix
	double **MIabs_med, **MIabs_uni_med, **MIabs_mult_med;// equal time correlation matrix

	double **TCS;						// matrix for Total  Correlation Strength
	double **CCS, **CCSmult;			// matrix for Cross- Correlation Strength w.r.t. uni and multivariate surrogates
	double **RCS, **RCSmult;			// matrix for Random Correlation Strength w.r.t. uni and multivariate surrogates

	double **TMI;						// matrix for Total  Mutual Information Strength
	double **CMI, **CMImult;			// matrix for Cross- Mutual Information Strength w.r.t. uni and multivariate surrogates
	double **RMI, **RMImult;			// matrix for Random Mutual Information Strength w.r.t. uni and multivariate surrogates

	int U;								// U statistic of U test
	double rho;							// Herrenstein's statistic for the overlap
	double p;							// significance

	double **CCS_p, **CCSmult_p;		// vector for sorting significances of matrix elements
	double **CMI_p, **CMImult_p;		// vector for sorting significances of matrix elements

// prototypes for Numerical Recipes routines

	unsigned char *cvector (long, long);
	int           *ivector (long, long);
	double         *vector (long, long);
	int          **imatrix (long, long, long, long);
	double        **matrix (long, long, long, long);
	double     ***f3tensor (long, long, long, long, long, long);

	void      free_cvector (unsigned char *, long, long);
	void      free_ivector      (int *, long, long);
	void       free_vector   (double *, long, long);
	void       free_matrix  (double **, long, long, long, long);
	void      free_imatrix     (int **, long, long, long, long);
	void    free_f3tensor  (double ***, long, long, long, long, long, long);

	void   sort (int, double[]);

// other prototypes

	void build_C        (int, int,      double, double **, double **);					// build C-Matrix
	void build_MI       (int, int, int, double, double **, double **, long);			// build MI matrix

	void Utest          (int, double[], int, double[], int *, double *, double *);		// Mann-Whitney-Wilcoxon U-test

	void IFS2 (int, int, double, int, double **, double ***, long, int, int, double, int **, double **);// generate Independent  Fourier based Surrogates (IAAFT)
	void MFS2 (int, int, double, int, double **, double ***, long, int, int, double, int *,  double **);// generate Multivariate Fourier based Surrogates (IAAFT)

	void Bonferroni_symmatr (int, double **, double **, double, int);					// Bonferroni correction of matrix elements
	void HolmBonf_symmatr   (int, double **, double **, double, int);					// family-wise error correction of matrix elements
	void BenjHoch_symmatr   (int, double **, double **, double, int);					// false discovery rate correction of matrix elements
	void uncorr_symmatr     (int, double **, double **, double, int);					// no multiple comparison correction

	printf ("  DIRECTORIES ../data and ../inputs MUST be provided!\n");
	printf ("  - ../data/series_raw.dat      is the multivariate time series to be analyzed\n");
	printf ("  - ../data                     all output files are written here\n");
	printf ("  - ../inputs/series_raw.inp    contains information about the time series\n");
	printf ("  - ../inputs/XMImatr_v2.inp    sets default parameters for the analysis\n");
	printf ("  EXIT with ctrl-C if files are inexistent or empty\n");
	printf ("\n");

// read input file for time series to be taken into account in multivariate analysis

	printf("\n");
    printf("***                 properties of TIME SERIES                     ***\n");
    printf("\n");

	fin_series = fopen("../inputs/series_raw_plosone.inp", "r");

	fscanf (fin_series, "%i",  &sw_time_in);
	if (sw_time_in!=0)
		printf ("  time given in column                                  %6i\n", sw_time_in);

	fscanf (fin_series, "%lf", &f_samp);
	printf ("  sampling frequency of recording                         %8.4f\n", f_samp);
	fscanf (fin_series, "%i",  &columns);
	printf ("  number of columns                                         %6i\n", columns);

	sel = ivector (1,columns);
	printf ("  columns included in analysis:  ");
	M = 0;
	for (ii=1;ii<=columns;ii++) {
		fscanf (fin_series, "%i", &sel[ii]);
		if (sel[ii]==1) {
			M++;
			printf(" %4i", ii);
		}
	}
	fclose(fin_series);

// read input file and set parameter values for calculation to default values

	printf("\n\n");
    printf("*** program XMImatr_v2_plosone analyzes multivariate time series  ***\n");
    printf("\n");

    printf("    PARAMETERS are set to the following values:\n\n");

	fin_calc = fopen("../inputs/XMImatr_v2_plosone.inp", "r");

	fscanf (fin_calc, "%i",  &sw_slopes);
	if (sw_slopes==1)
		printf ("  first temporal derivatives are analyzed instead of original signals!\n");

	fscanf (fin_calc, "%i",  &neighbors);
	printf ("  XCS matrices are calculated on basis of equal-time cross-correlation\n");
	if (neighbors>0) {
		printf ("  XMI matrices are calculated on basis of mutual information\n");
		printf ("  number of nearest neighbors            k =                %6i\n", neighbors);
	}

	fscanf (fin_calc, "%i",  &T_surr);
	printf ("  length of moving window for surrogates T_surr =           %6i\n", T_surr);
	fscanf (fin_calc, "%i",  &T);
	printf ("  length of moving window for matrices   T =                %6i\n", T);
	fscanf (fin_calc, "%i",  &N_ens);
	printf ("  ensemble size for statistics           N_ens =            %6i\n", N_ens);
	fscanf (fin_calc, "%i",  &N_surr);
	printf ("  number of surrogates                   N_surr =           %6i\n", N_surr);
	fscanf (fin_calc, "%i",  &step_surr);
	printf ("  step width of moving window                               %6i\n", step_surr);
	printf ("\n");

	fscanf (fin_calc, "%lf", &tstart);
	if (tstart<1.0*T_surr/f_samp)
		tstart = 1.0*T_surr/f_samp;
	printf ("  analysis           starts     at second                 %8.2f\n", tstart);
	fscanf (fin_calc, "%lf", &tstop);
	printf ("  analysis           terminates at second                 %8.2f\n", tstop);
	fscanf (fin_calc, "%lf", &siglev);
	printf ("  significance level for statistics                       %8.2f\n", siglev);
	fscanf (fin_calc, "%i", &sw_corr);
	if (sw_corr==0)
		printf ("  CCS/CMI: uncorrected matrix elements (expected to produce many false positives!)\n");
	if (sw_corr==1)
		printf ("  CCS/CMI: Benjamini-Hochberg / false discovery rate correction\n");
	if (sw_corr==2)
		printf ("  CCS/CMI: Holm-Bonferroni / family-wise error rate correction\n");
	if (sw_corr==3)
		printf ("  CCS/CMI: Bonferroni correction (expected to produce many false negatives!)\n");
	printf ("\n");

    printf("\n  other SETTINGS:\n\n");

	fscanf (fin_calc, "%i",  &sw_time_out);
	if (sw_time_out!=0)
		printf ("  time is written to output file in columns               1 ... %2i\n", sw_time_out);
	printf ("\n");

	fclose (fin_calc);

// name output files

	if (sw_slopes==0)
		strcat (name_fout_matrix_TCS, "TCS_matrix.dat");
	else
		strcat (name_fout_matrix_TCS, "STCS_matrix.dat");
	printf ("  output file for univariate   TCS matrix    %s\n", name_fout_matrix_TCS);

	if (sw_slopes==0)
		strcat (name_fout_matrix_CCS, "uniCCSpos_matrix.dat");
	else
		strcat (name_fout_matrix_CCS, "uniSCCSpos_matrix.dat");
	printf ("  output file for univariate   CCS matrix    %s\n", name_fout_matrix_CCS);

	if (sw_slopes==0)
		strcat (name_fout_matrix_RCS, "uniRCS_matrix.dat");
	else
		strcat (name_fout_matrix_RCS, "uniSRCS_matrix.dat");
	printf ("  output file for univariate   RCS matrix    %s\n", name_fout_matrix_RCS);

	if (neighbors>0) {
		if (sw_slopes==0)
			strcat (name_fout_matrix_multCCS, "multCCSpos_matrix.dat");
		else
			strcat (name_fout_matrix_multCCS, "multSCCSpos_matrix.dat");
		printf ("  output file for multivariate CCS matrix    %s\n", name_fout_matrix_multCCS);

		if (sw_slopes==0)
			strcat (name_fout_matrix_multRCS, "multRCS_matrix.dat");
		else
			strcat (name_fout_matrix_multRCS, "multSRCS_matrix.dat");
		printf ("  output file for multivariate RCS matrix    %s\n", name_fout_matrix_multRCS);

		printf ("\n");

		if (sw_slopes==0)
			strcat (name_fout_matrix_TMI, "TMI_matrix.dat");
		else
			strcat (name_fout_matrix_TMI, "STMI_matrix.dat");
		printf ("  output file for univariate   TMI matrix    %s\n", name_fout_matrix_TMI);

		if (sw_slopes==0)
			strcat (name_fout_matrix_CMI, "uniCMIpos_matrix.dat");
		else
			strcat (name_fout_matrix_CMI, "uniSCMIpos_matrix.dat");
		printf ("  output file for univariate   CMI matrix    %s\n", name_fout_matrix_CMI);

		if (sw_slopes==0)
			strcat (name_fout_matrix_RMI, "uniRMI_matrix.dat");
		else
			strcat (name_fout_matrix_RMI, "uniSRMI_matrix.dat");
		printf ("  output file for univariate   RMI matrix    %s\n", name_fout_matrix_RMI);

		if (sw_slopes==0)
			strcat (name_fout_matrix_multCMI, "multCMIpos_matrix.dat");
		else
			strcat (name_fout_matrix_multCMI, "multSCMIpos_matrix.dat");
		printf ("  output file for multivariate CMI matrix    %s\n", name_fout_matrix_multCMI);

		if (sw_slopes==0)
			strcat (name_fout_matrix_multRMI, "multRMI_matrix.dat");
		else
			strcat (name_fout_matrix_multRMI, "multSRMI_matrix.dat");
		printf ("  output file for multivariate RMI matrix    %s\n", name_fout_matrix_multRMI);
	}

// open output files

	fout_matrix_TCS = fopen(name_fout_matrix_TCS, "w");
	fprintf(fout_matrix_TCS, "* output file of program XMImatr_v2_plosone.c\n");
	fprintf(fout_matrix_TCS, "*\n");
	fprintf(fout_matrix_TCS, "* elements of TCS matrix of multivariate time series\n");
	if (sw_slopes==1)
		fprintf(fout_matrix_TCS, "* first temporal derivatives are analyzed instead of original signals!\n");
	fprintf(fout_matrix_TCS, "* window lengths: for surrogates T_surr = %6i, for interrelation matrices T = %6i\n", T_surr, T);
	fprintf(fout_matrix_TCS, "* significance level for statistical tests alpha = %8.4f\n", siglev);
	fprintf(fout_matrix_TCS, "*\n");

	fout_matrix_CCS = fopen(name_fout_matrix_CCS, "w");
	fprintf(fout_matrix_CCS, "* output file of program XMImatr_v2_plosone.c\n");
	fprintf(fout_matrix_CCS, "*\n");
	fprintf(fout_matrix_CCS, "* elements of CCS matrix of multivariate time series\n");
	if (sw_slopes==1)
		fprintf(fout_matrix_CCS, "* first temporal derivatives are analyzed instead of original signals!\n");
	fprintf(fout_matrix_CCS, "* window lengths: for surrogates T_surr = %6i, for interrelation matrices T = %6i\n", T_surr, T);
	fprintf(fout_matrix_CCS, "* significance level for statistical tests alpha = %8.4f\n", siglev);
	if (sw_corr==0)
		fprintf (fout_matrix_CCS, "* uncorrected matrix elements (expected to produce many false positives!)\n");
	if (sw_corr==1)
		fprintf (fout_matrix_CCS, "* Benjamini-Hochberg / false discovery rate correction\n");
	if (sw_corr==2)
		fprintf (fout_matrix_CCS, "* Holm-Bonferroni / family-wise error rate correction\n");
	if (sw_corr==3)
		fprintf (fout_matrix_CCS, "* Bonferroni correction with factor %6i (expected to produce many false negatives!)\n", M*(M-1)/2);
	fprintf(fout_matrix_CCS, "*\n");

	fout_matrix_RCS = fopen(name_fout_matrix_RCS, "w");
	fprintf(fout_matrix_RCS, "* output file of program XMImatr_v2_plosone.c\n");
	fprintf(fout_matrix_RCS, "*\n");
	fprintf(fout_matrix_RCS, "* elements of RCS matrix of multivariate time series\n");
	if (sw_slopes==1)
		fprintf(fout_matrix_RCS, "* first temporal derivatives are analyzed instead of original signals!\n");
	fprintf(fout_matrix_RCS, "* window lengths: for surrogates T_surr = %6i, for interrelation matrices T = %6i\n", T_surr, T);
	fprintf(fout_matrix_RCS, "* significance level for statistical tests alpha = %8.4f\n", siglev);
	fprintf(fout_matrix_RCS, "*\n");

	if (neighbors>0) {
		fout_matrix_multCCS = fopen(name_fout_matrix_multCCS, "w");
		fprintf(fout_matrix_multCCS, "* output file of program XMImatr_v2_plosone.c\n");
		fprintf(fout_matrix_multCCS, "*\n");
		fprintf(fout_matrix_multCCS, "* elements of CCS matrix of multivariate time series\n");
		if (sw_slopes==1)
			fprintf(fout_matrix_multCCS, "* first temporal derivatives are analyzed instead of original signals!\n");
		fprintf(fout_matrix_multCCS, "* window lengths: for surrogates T_surr = %6i, for interrelation matrices T = %6i\n", T_surr, T);
		fprintf(fout_matrix_multCCS, "* significance level for statistical tests alpha = %8.4f\n", siglev);
		if (sw_corr==0)
			fprintf (fout_matrix_multCCS, "* uncorrected matrix elements (expected to produce many false positives!)\n");
		if (sw_corr==1)
			fprintf (fout_matrix_multCCS, "* Benjamini-Hochberg / false discovery rate correction\n");
		if (sw_corr==2)
			fprintf (fout_matrix_multCCS, "* Holm-Bonferroni / family-wise error rate correction\n");
		if (sw_corr==3)
			fprintf (fout_matrix_multCCS, "* Bonferroni correction with factor %6i (expected to produce many false negatives!)\n", M*(M-1)/2);
		fprintf(fout_matrix_multCCS, "*\n");

		fout_matrix_multRCS = fopen(name_fout_matrix_multRCS, "w");
		fprintf(fout_matrix_multRCS, "* output file of program XMImatr_v2_plosone.c\n");
		fprintf(fout_matrix_multRCS, "*\n");
		fprintf(fout_matrix_multRCS, "* elements of RCS matrix of multivariate time series\n");
		if (sw_slopes==1)
			fprintf(fout_matrix_multRCS, "* first temporal derivatives are analyzed instead of original signals!\n");
		fprintf(fout_matrix_multRCS, "* window lengths: for surrogates T_surr = %6i, for interrelation matrices T = %6i\n", T_surr, T);
		fprintf(fout_matrix_multRCS, "* significance level for statistical tests alpha = %8.4f\n", siglev);
		fprintf(fout_matrix_multRCS, "*\n");

		fout_matrix_TMI = fopen(name_fout_matrix_TMI, "w");
		fprintf(fout_matrix_TMI, "* output file of program XMImatr_v2_plosone.c\n");
		fprintf(fout_matrix_TMI, "*\n");
		fprintf(fout_matrix_TMI, "* elements of TMI matrix of multivariate time series, number of nearest neighbors k = %4i\n", neighbors);
		if (sw_slopes==1)
			fprintf(fout_matrix_TMI, "* first temporal derivatives are analyzed instead of original signals!\n");
		fprintf(fout_matrix_TMI, "* window lengths: for surrogates T_surr = %6i, for interrelation matrices T = %6i\n", T_surr, T);
		fprintf(fout_matrix_TMI, "* significance level for statistical tests alpha = %8.4f\n", siglev);
		fprintf(fout_matrix_TMI, "*\n");

		fout_matrix_CMI = fopen(name_fout_matrix_CMI, "w");
		fprintf(fout_matrix_CMI, "* output file of program XMImatr_v2_plosone.c\n");
		fprintf(fout_matrix_CMI, "*\n");
		fprintf(fout_matrix_CMI, "* elements of CMI matrix of multivariate time series, number of nearest neighbors k = %4i\n", neighbors);
		if (sw_slopes==1)
			fprintf(fout_matrix_CMI, "* first temporal derivatives are analyzed instead of original signals!\n");
		fprintf(fout_matrix_CMI, "* window lengths: for surrogates T_surr = %6i, for interrelation matrices T = %6i\n", T_surr, T);
		fprintf(fout_matrix_CMI, "* significance level for statistical tests alpha = %8.4f\n", siglev);
		if (sw_corr==0)
			fprintf (fout_matrix_CMI, "* uncorrected matrix elements (expected to produce many false positives!)\n");
		if (sw_corr==1)
			fprintf (fout_matrix_CMI, "* Benjamini-Hochberg / false discovery rate correction\n");
		if (sw_corr==2)
			fprintf (fout_matrix_CMI, "* Holm-Bonferroni / family-wise error rate correction\n");
		if (sw_corr==3)
			fprintf (fout_matrix_CMI, "* Bonferroni correction with factor %6i (expected to produce many false negatives!)\n", M*(M-1)/2);
		fprintf(fout_matrix_CMI, "*\n");

		fout_matrix_RMI = fopen(name_fout_matrix_RMI, "w");
		fprintf(fout_matrix_RMI, "* output file of program XMImatr_v2_plosone.c\n");
		fprintf(fout_matrix_RMI, "*\n");
		fprintf(fout_matrix_RMI, "* elements of RMI matrix of multivariate time series, number of nearest neighbors k = %4i\n", neighbors);
		if (sw_slopes==1)
			fprintf(fout_matrix_RMI, "* first temporal derivatives are analyzed instead of original signals!\n");
		fprintf(fout_matrix_RMI, "* window lengths: for surrogates T_surr = %6i, for interrelation matrices T = %6i\n", T_surr, T);
		fprintf(fout_matrix_RMI, "* significance level for statistical tests alpha = %8.4f\n", siglev);
		fprintf(fout_matrix_RMI, "*\n");

		fout_matrix_multCMI = fopen(name_fout_matrix_multCMI, "w");
		fprintf(fout_matrix_multCMI, "* output file of program XMImatr_v2_plosone.c\n");
		fprintf(fout_matrix_multCMI, "*\n");
		fprintf(fout_matrix_multCMI, "* elements of CMI matrix of multivariate time series, number of nearest neighbors k = %4i\n", neighbors);
		if (sw_slopes==1)
			fprintf(fout_matrix_multCMI, "* first temporal derivatives are analyzed instead of original signals!\n");
		fprintf(fout_matrix_multCMI, "* window lengths: for surrogates T_surr = %6i, for interrelation matrices T = %6i\n", T_surr, T);
		fprintf(fout_matrix_multCMI, "* significance level for statistical tests alpha = %8.4f\n", siglev);
		if (sw_corr==0)
			fprintf (fout_matrix_multCMI, "* uncorrected matrix elements (expected to produce many false positives!)\n");
		if (sw_corr==1)
			fprintf (fout_matrix_multCMI, "* Benjamini-Hochberg / false discovery rate correction\n");
		if (sw_corr==2)
			fprintf (fout_matrix_multCMI, "* Holm-Bonferroni / family-wise error rate correction\n");
		if (sw_corr==3)
			fprintf (fout_matrix_multCMI, "* Bonferroni correction with factor %6i (expected to produce many false negatives!)\n", M*(M-1)/2);
		fprintf(fout_matrix_multCMI, "*\n");

		fout_matrix_multRMI = fopen(name_fout_matrix_multRMI, "w");
		fprintf(fout_matrix_multRMI, "* output file of program XMImatr_v2_plosone.c\n");
		fprintf(fout_matrix_multRMI, "*\n");
		fprintf(fout_matrix_multRMI, "* elements of RMI matrix of multivariate time series, number of nearest neighbors k = %4i\n", neighbors);
		if (sw_slopes==1)
			fprintf(fout_matrix_multRMI, "* first temporal derivatives are analyzed instead of original signals!\n");
		fprintf(fout_matrix_multRMI, "* window lengths: for surrogates T_surr = %6i, for interrelation matrices T = %6i\n", T_surr, T);
		fprintf(fout_matrix_multRMI, "* significance level for statistical tests alpha = %8.4f\n", siglev);
		fprintf(fout_matrix_multRMI, "*\n");
	}

// allocation of arrays

	Xl             = vector            (1,M);
	Xl_old         = vector            (1,M);

	X_read         = matrix (0,T_surr-1,1,M);
	X_surr         = matrix (0,T_surr-1,1,M);

	X              = matrix      (0,T-1,1,M);
	Xtilde         = matrix      (0,T-1,1,M);

	X_IFS          =  f3tensor (0,T_surr-1,1,M,1,N_surr);
	X_IFSsing      =  matrix   (0,T_surr-1,1,M);
	C_IFS          =  matrix              (1,M,1,M);
	MI_IFS         =  matrix              (1,M,1,M);
	iter_uni       = imatrix              (1,M,1,N_surr);
	accur_uni      =  matrix              (1,M,1,N_surr);

	X_MFS          =  f3tensor (0,T_surr-1,1,M,1,N_surr);
	X_MFSsing      =  matrix   (0,T_surr-1,1,M);
	C_MFS          =  matrix              (1,M,1,M);
	MI_MFS         =  matrix              (1,M,1,M);
	iter_mult      = ivector                  (1,N_surr);
	accur_mult     =  matrix              (1,M,1,N_surr);

	C              = matrix        (1,M,1,M);
	all_C          = f3tensor      (1,M,1,M,1,N_ens);
	all_Cabs       = f3tensor      (1,M,1,M,1,N_ens);
	sort_C         = vector                (1,N_ens);
	sort_Cabs      = vector                (1,N_ens);
	C_med          = matrix        (1,M,1,M);
	Cabs_med       = matrix        (1,M,1,M);

	vec_med        = vector                (1,N_ens);

	C_surr         = matrix        (1,M,1,M);
	all_C_uni      = f3tensor      (1,M,1,M,1,N_ens*N_surr);
	all_Cabs_uni   = f3tensor      (1,M,1,M,1,N_ens*N_surr);
	all_C_mult     = f3tensor      (1,M,1,M,1,N_ens*N_surr);
	all_Cabs_mult  = f3tensor      (1,M,1,M,1,N_ens*N_surr);
	sort_C_surr    = vector                (1,N_ens*N_surr);
	sort_Cabs_uni  = vector                (1,N_ens*N_surr);
	sort_Cabs_mult = vector                (1,N_ens*N_surr);
	C_uni_med      = matrix        (1,M,1,M);
	Cabs_uni_med   = matrix        (1,M,1,M);
	C_mult_med     = matrix        (1,M,1,M);
	Cabs_mult_med  = matrix        (1,M,1,M);

	TCS            = matrix        (1,M,1,M);
	CCS            = matrix        (1,M,1,M);
	RCS            = matrix        (1,M,1,M);
	CCSmult        = matrix        (1,M,1,M);
	RCSmult        = matrix        (1,M,1,M);

	CCS_p          = matrix        (1,M,1,M);
	CCSmult_p      = matrix        (1,M,1,M);

	if (neighbors>0) {
		MI             = matrix        (1,M,1,M);
		all_MI         = f3tensor      (1,M,1,M,1,N_ens);
		all_MIabs      = f3tensor      (1,M,1,M,1,N_ens);
		sort_MI        = vector                (1,N_ens);
		sort_MIabs     = vector                (1,N_ens);
		MI_med         = matrix        (1,M,1,M);
		MIabs_med      = matrix        (1,M,1,M);

		MI_surr        = matrix        (1,M,1,M);
		all_MI_uni     = f3tensor      (1,M,1,M,1,N_ens*N_surr);
		all_MIabs_uni  = f3tensor      (1,M,1,M,1,N_ens*N_surr);
		all_MI_mult    = f3tensor      (1,M,1,M,1,N_ens*N_surr);
		all_MIabs_mult = f3tensor      (1,M,1,M,1,N_ens*N_surr);
		sort_MI_surr   = vector                (1,N_ens*N_surr);
		sort_MIabs_uni = vector                (1,N_ens*N_surr);
		sort_MIabs_mult= vector                (1,N_ens*N_surr);
		MI_uni_med     = matrix        (1,M,1,M);
		MIabs_uni_med  = matrix        (1,M,1,M);
		MI_mult_med    = matrix        (1,M,1,M);
		MIabs_mult_med = matrix        (1,M,1,M);

		TMI            = matrix        (1,M,1,M);
		CMI            = matrix        (1,M,1,M);
		RMI            = matrix        (1,M,1,M);
		CMImult        = matrix        (1,M,1,M);
		RMImult        = matrix        (1,M,1,M);

		CMI_p          = matrix        (1,M,1,M);
		CMImult_p      = matrix        (1,M,1,M);
	}

// *** estimation of interrelation matrices from multivariate time series ***

	printf("\n");
    printf("*** ANALYZING time series ...                                     ***\n");
	printf("\n");
	printf("    progress:\n");

	fin_series = fopen("../data/series_raw.dat", "r");

	step = T - (int)((N_ens*T-T_surr)/(N_ens-1)+1);

	(void) time (&t1);
	seed = (long) t1; seed = -seed;

	line = cvector (0,columns*30);

	l = 0;
	while (fgets(line,columns*30,fin_series)!=NULL) {
		if (line[0]!='*') {
			l++;

// *** outer moving window for surrogate generation ***

// continuously fill matrix X_read covering moving time window of length T_surr

			i = 0;
			for (ii=1;ii<=columns;ii++) {
				fscanf(fin_series, "%lf", &readdata);
				if (ii==sw_time_in)
					t = readdata;

				if (sel[ii]==1) {
					i++;
					if (sw_slopes==0)
						X_read[l%T_surr][i] = readdata;
					else
						X_read[l%T_surr][i] = readdata - X_read[(l-1)%T_surr][i];
				}
			}

			if (sw_time_in==0)
//			if ((sw_time_in==0) || (M==2))			// why necessary ???
				t = 1.0*l/f_samp;

			if ((t==(int)(t)) && (t>=tstart))
				printf("  %6i of %6i seconds done\n", (int)(t-tstart), (int)(tstop-tstart));
			if (t+T_surr/f_samp<tstart)
				continue;
			if (t>tstop)
				break;

// for multiples of step_surr

			if ((t>=tstart) && (l%step_surr==0)) {

// fill matrix X_surr used for surrogate generation (outer window)

				for (k=1;k<=T_surr;k++) {
					for (i=1;i<=M;i++)
						X_surr[k-1][i] = X_read[(l+k)%T_surr][i];
				}

// generate ensemble of univariate surrogates

				IFS2 (M,T_surr,f_samp, N_surr,X_surr,X_IFS, seed-l,1, 500,0.0001, iter_uni,accur_uni);

// write example time series to output files
/*
				fout_TEST = fopen("../data/test_EEG.dat", "w");
				for (k=1;k<=T_surr;k++) {
					for (i=1;i<=M;i++)
						fprintf(fout_TEST, "  %12.4e", X_surr[k-1][i]);
					fprintf(fout_TEST, "\n");
				}
				fclose(fout_TEST);

				fout_TEST = fopen("../data/test_IFS1.dat", "w");
				for (k=1;k<=T_surr;k++) {
					for (i=1;i<=M;i++)
						fprintf(fout_TEST, "  %12.4e", X_IFS[k-1][i][1]);
					fprintf(fout_TEST, "\n");
				}
				fclose(fout_TEST);

				fout_TEST = fopen("../data/test_IFS2.dat", "w");
				for (k=1;k<=T_surr;k++) {
					for (i=1;i<=M;i++)
						fprintf(fout_TEST, "  %12.4e", X_IFS[k-1][i][N_surr]);
					fprintf(fout_TEST, "\n");
				}
				fclose(fout_TEST);

				for (n=1;n<=N_ens;n++)
					printf("  %6i", iter_uni[2][n]);
				printf("\n");

				for (n=1;n<=N_ens;n++)
					printf("  %12.4e", accur_uni[2][n]);
				printf("\n");
*/
// optionally generate ensemble of multivariate surrogates

				if (neighbors>0)
					MFS2 (M,T_surr,f_samp, N_surr,X_surr,X_MFS, seed-l,1, 1000,0.01, iter_mult,accur_mult);

// write example time series to output files
/*
				fout_TEST = fopen("../data/test_MFS1.dat", "w");
				for (k=1;k<=T_surr;k++) {
					for (i=1;i<=M;i++)
						fprintf(fout_TEST, "  %12.4e", X_MFS[k-1][i][1]);
					fprintf(fout_TEST, "\n");
				}
				fclose(fout_TEST);

				fout_TEST = fopen("../data/test_MFS2.dat", "w");
				for (k=1;k<=T_surr;k++) {
					for (i=1;i<=M;i++)
						fprintf(fout_TEST, "  %12.4e", X_MFS[k-1][i][N_surr]);
					fprintf(fout_TEST, "\n");
				}
				fclose(fout_TEST);

				for (n=1;n<=N_ens;n++)
					printf("  %6i", iter_mult[n]);
				printf("\n");

				for (n=1;n<=N_ens;n++)
					printf("  %12.4e", accur_mult[2][n]);
				printf("\n");
*/
// *** inner moving window ***

				for (n=0;n<N_ens;n++) {

// *** calculation for original data ***

// normalize windows of length T and fill data matrix

					for (i=1;i<=M;i++) {
						ave = 0.0;
						dev = 0.0;
						for (k=0;k<T;k++) {
							ave = ave + X_surr[n*step+k][i];
							dev = dev + X_surr[n*step+k][i]*X_surr[n*step+k][i];
						}
						ave =       ave/T;
						dev = sqrt (dev/T - ave*ave);

						for (k=0;k<T;k++)
							Xtilde[k][i] = (X_surr[n*step+k][i] - ave) / dev;
					}

// fill interrelation matrices

					build_C (M,T,1.0,Xtilde,C);
					if (neighbors>0)
						build_MI (M,T,neighbors,1.0,Xtilde,MI,seed-l);

					for (i=1;i<=M;i++) {
						for (j=1;j<=M;j++) {
							all_C   [i][j][n] =      C[i][j];
							all_Cabs[i][j][n] = fabs(C[i][j]);

							if (neighbors>0) {
								all_MI   [i][j][n] =      MI[i][j];
								all_MIabs[i][j][n] = fabs(MI[i][j]);
							}
						}
					}

// *** calculation for univariate surrogate data ***

					for (nn=1;nn<=N_surr;nn++) {

// normalize windows of length T and fill data matrix

						for (i=1;i<=M;i++) {
							ave = 0.0;
							dev = 0.0;
							for (k=0;k<T;k++) {
								ave = ave + X_IFS[n*step+k][i][nn];
								dev = dev + X_IFS[n*step+k][i][nn]*X_IFS[n*step+k][i][nn];
							}
							ave =       ave/T;
							dev = sqrt (dev/T - ave*ave);

							for (k=0;k<T;k++)
								Xtilde[k][i] = (X_IFS[n*step+k][i][nn] - ave) / dev;
						}

// fill interrelation matrices

						build_C (M,T,1.0,Xtilde,C_surr);
						if (neighbors>0)
							build_MI (M,T,neighbors,1.0,Xtilde,MI_surr,seed-l);

						for (i=1;i<=M;i++) {
							for (j=1;j<=M;j++) {
								all_C_uni   [i][j][n*N_surr+nn] =      C_surr[i][j];
								all_Cabs_uni[i][j][n*N_surr+nn] = fabs(C_surr[i][j]);

								if (neighbors>0) {
									all_MI_uni   [i][j][n*N_surr+nn] =      MI_surr[i][j];
									all_MIabs_uni[i][j][n*N_surr+nn] = fabs(MI_surr[i][j]);
								}
							}
						}
					}

// *** calculation for multivariate surrogate data ***

					if (neighbors>0) {
						for (nn=1;nn<=N_surr;nn++) {

// normalize windows of length T and fill data matrix

							for (i=1;i<=M;i++) {
								ave = 0.0;
								dev = 0.0;
								for (k=0;k<T;k++) {
									ave = ave + X_MFS[n*step+k][i][nn];
									dev = dev + X_MFS[n*step+k][i][nn]*X_MFS[n*step+k][i][nn];
								}
								ave =       ave/T;
								dev = sqrt (dev/T - ave*ave);

								for (k=0;k<T;k++)
									Xtilde[k][i] = (X_MFS[n*step+k][i][nn] - ave) / dev;
							}

// fill interrelation matrices

							build_C (M,T,1.0,Xtilde,C_surr);
							build_MI (M,T,neighbors,1.0,Xtilde,MI_surr,seed-l);

							for (i=1;i<=M;i++) {
								for (j=1;j<=M;j++) {
									all_C_mult   [i][j][n*N_surr+nn] =      C_surr[i][j];
									all_Cabs_mult[i][j][n*N_surr+nn] = fabs(C_surr[i][j]);

									if (neighbors>0) {
										all_MI_mult   [i][j][n*N_surr+nn] =      MI_surr[i][j];
										all_MIabs_mult[i][j][n*N_surr+nn] = fabs(MI_surr[i][j]);
									}
								}
							}
						}
					}
				}

// *** calculation of matrices TCS, CCS and RCS ***
// references:
// Rummel et al., J. Neurosci. Meth. 191, 94-100 (2010)
// Rummel et al., Phys. Rev. E83, 066215 (2011)

				for (i=1;i<=M;i++) {
					TCS    [i][i] = 1.0;
					RCS    [i][i] = 1.0;
					CCS    [i][i] = 1.0;
					RCSmult[i][i] = 1.0;
					CCSmult[i][i] = 1.0;

					for (j=i+1;j<=M;j++) {

// medians of original matrix elements: TCS matrix

						for (n=1;n<=N_ens;n++) {
							sort_C   [n] = all_C   [i][j][n];
							sort_Cabs[n] = all_Cabs[i][j][n];
						}
						sort (N_ens, sort_C);
						sort (N_ens, sort_Cabs);
						if (0.5*N_ens!=N_ens/2) {
							C_med   [i][j] = sort_C   [N_ens/2+1];
							Cabs_med[i][j] = sort_Cabs[N_ens/2+1];
						}
						else {
							C_med   [i][j] = 0.5*(sort_C   [N_ens/2]+sort_C   [N_ens/2+1]);
							Cabs_med[i][j] = 0.5*(sort_Cabs[N_ens/2]+sort_Cabs[N_ens/2+1]);
						}

						if (C_med[i][j]>0.0)
							TCS[i][j] =  Cabs_med[i][j];
						else
							TCS[i][j] = -Cabs_med[i][j];

						TCS[j][i] = TCS[i][j];

// medians of phase randomized surrogate matrix elements: RCS matrices

						for (n=1;n<=N_ens*N_surr;n++)
							sort_Cabs_uni[n] = all_Cabs_uni[i][j][n];
						sort (N_ens*N_surr, sort_Cabs_uni);
						if (0.5*N_ens*N_surr!=N_ens*N_surr/2)
							Cabs_uni_med[i][j] = sort_Cabs_uni[N_ens*N_surr/2+1];
						else
							Cabs_uni_med[i][j] = 0.5*(sort_Cabs_uni[N_ens*N_surr/2]+sort_Cabs_uni[N_ens*N_surr/2+1]);

						if (C_med[i][j]>0.0)
							RCS[i][j] =  Cabs_uni_med[i][j];
						else
							RCS[i][j] = -Cabs_uni_med[i][j];

						RCS[j][i] = RCS[i][j];

						if (neighbors>0) {
							for (n=1;n<=N_ens*N_surr;n++)
								sort_Cabs_mult[n] = all_Cabs_mult[i][j][n];
							sort (N_ens*N_surr, sort_Cabs_mult);
							if (0.5*N_ens*N_surr!=N_ens*N_surr/2)
								Cabs_mult_med[i][j] = sort_Cabs_mult[N_ens*N_surr/2+1];
							else
								Cabs_mult_med[i][j] = 0.5*(sort_Cabs_mult[N_ens*N_surr/2]+sort_Cabs_mult[N_ens*N_surr/2+1]);

							if (C_med[i][j]>0.0)
								RCSmult[i][j] =  Cabs_mult_med[i][j];
							else
								RCSmult[i][j] = -Cabs_mult_med[i][j];

							RCSmult[j][i] = RCSmult[i][j];
						}

// define CCS matrices
// optional correction for multiple tests is performed afterwards

						Utest (N_ens,sort_Cabs, N_ens*N_surr,sort_Cabs_uni, &U,&rho,&p);
						CCS_p[i][j] = p;

						if (Cabs_med[i][j]>Cabs_uni_med[i][j]) {
							if (C_med[i][j]>0.0)
								CCS[i][j] =  (Cabs_med[i][j] - Cabs_uni_med[i][j]) / (1.0 - Cabs_uni_med[i][j]);
							else
								CCS[i][j] = -(Cabs_med[i][j] - Cabs_uni_med[i][j]) / (1.0 - Cabs_uni_med[i][j]);
						}
						else
							CCS[i][j] = 0.0;

						if (neighbors>0) {
							Utest (N_ens,sort_Cabs, N_ens*N_surr,sort_Cabs_mult, &U,&rho,&p);
							CCSmult_p[i][j] = p;

							if (Cabs_med[i][j]>Cabs_mult_med[i][j]) {
								if (C_med[i][j]>0.0)
									CCSmult[i][j] =  (Cabs_med[i][j] - Cabs_mult_med[i][j]) / (1.0 - Cabs_mult_med[i][j]);
								else
									CCSmult[i][j] = -(Cabs_med[i][j] - Cabs_mult_med[i][j]) / (1.0 - Cabs_mult_med[i][j]);
							}
							else
								CCSmult[i][j] = 0.0;
						}
					}
				}

// *** calculation of matrices TMI, CMI and RMI ***
// reference:
// Rummel et al., Phys. Rev. E83, 066215 (2011)

				if (neighbors>0) {
					for (i=1;i<=M;i++) {
						TMI    [i][i] = 1.0;
						RMI    [i][i] = 1.0;
						CMI    [i][i] = 1.0;
						RMImult[i][i] = 1.0;
						CMImult[i][i] = 1.0;

						for (j=i+1;j<=M;j++) {

// medians of original matrix elements: TMI matrix

							for (n=1;n<=N_ens;n++) {
								sort_MI   [n] = all_MI   [i][j][n];
								sort_MIabs[n] = all_MIabs[i][j][n];
							}
							sort (N_ens, sort_MI);
							sort (N_ens, sort_MIabs);
							if (0.5*N_ens!=N_ens/2) {
								MI_med   [i][j] = sort_MI   [N_ens/2+1];
								MIabs_med[i][j] = sort_MIabs[N_ens/2+1];
							}
							else {
								MI_med   [i][j] = 0.5*(sort_MI   [N_ens/2]+sort_MI   [N_ens/2+1]);
								MIabs_med[i][j] = 0.5*(sort_MIabs[N_ens/2]+sort_MIabs[N_ens/2+1]);
							}

							if (MI_med[i][j]>0.0)
								TMI[i][j] =  MIabs_med[i][j];
							else
								TMI[i][j] = -MIabs_med[i][j];

							TMI[j][i] = TMI[i][j];

// medians of phase randomized surrogate matrix elements: RMI matrices

							for (n=1;n<=N_ens*N_surr;n++)
								sort_MIabs_uni[n] = all_MIabs_uni[i][j][n];
							sort (N_ens*N_surr, sort_MIabs_uni);
							if (0.5*N_ens*N_surr!=N_ens*N_surr/2)
								MIabs_uni_med[i][j] = sort_MIabs_uni[N_ens*N_surr/2+1];
							else
								MIabs_uni_med[i][j] = 0.5*(sort_MIabs_uni[N_ens*N_surr/2]+sort_MIabs_uni[N_ens*N_surr/2+1]);

							if (MI_med[i][j]>0.0)
								RMI[i][j] =  MIabs_uni_med[i][j];
							else
								RMI[i][j] = -MIabs_uni_med[i][j];

							RMI[j][i] = RMI[i][j];

							for (n=1;n<=N_ens*N_surr;n++)
								sort_MIabs_mult[n] = all_MIabs_mult[i][j][n];
							sort (N_ens*N_surr, sort_MIabs_mult);
							if (0.5*N_ens*N_surr!=N_ens*N_surr/2)
								MIabs_mult_med[i][j] = sort_MIabs_mult[N_ens*N_surr/2+1];
							else
								MIabs_mult_med[i][j] = 0.5*(sort_MIabs_mult[N_ens*N_surr/2]+sort_MIabs_mult[N_ens*N_surr/2+1]);

							if (MI_med[i][j]>0.0)
								RMImult[i][j] =  MIabs_mult_med[i][j];
							else
								RMImult[i][j] = -MIabs_mult_med[i][j];

							RMImult[j][i] = RMImult[i][j];

// define CMI matrices
// optional correction for multiple tests is performed afterwards

							Utest (N_ens,sort_MIabs, N_ens*N_surr,sort_MIabs_uni, &U,&rho,&p);
							CMI_p[i][j] = p;

							if (MIabs_med[i][j]>MIabs_uni_med[i][j]) {
								if (MI_med[i][j]>0.0)
									CMI[i][j] =  (MIabs_med[i][j] - MIabs_uni_med[i][j]) / (1.0 - MIabs_uni_med[i][j]);
								else
									CMI[i][j] = -(MIabs_med[i][j] - MIabs_uni_med[i][j]) / (1.0 - MIabs_uni_med[i][j]);
							}
							else
								CMI[i][j] = 0.0;

							Utest (N_ens,sort_MIabs, N_ens*N_surr,sort_MIabs_mult, &U,&rho,&p);
							CMImult_p[i][j] = p;

							if (MIabs_med[i][j]>MIabs_mult_med[i][j]) {
								if (MI_med[i][j]>0.0)
									CMImult[i][j] =  (MIabs_med[i][j] - MIabs_mult_med[i][j]) / (1.0 - MIabs_mult_med[i][j]);
								else
									CMImult[i][j] = -(MIabs_med[i][j] - MIabs_mult_med[i][j]) / (1.0 - MIabs_mult_med[i][j]);
							}
							else
								CMImult[i][j] = 0.0;
						}
					}
				}

// *** delete insignificant matrix elements and optionally correct statistics for multiple comparisons

// no multiple comparison correction (expected to produce many false positives!)

				if (sw_corr==0) {
					uncorr_symmatr (M, CCS,     CCS_p,     siglev, 0);
					uncorr_symmatr (M, CCSmult, CCSmult_p, siglev, 0);

					if (neighbors>0) {
						uncorr_symmatr (M, CMI,     CMI_p,     siglev, 0);
						uncorr_symmatr (M, CMImult, CMImult_p, siglev, 0);
					}
				}

// Benjamini-Hochberg / false discovery rate correction

				if (sw_corr==1) {
					BenjHoch_symmatr (M, CCS,     CCS_p,     siglev, 0);
					BenjHoch_symmatr (M, CCSmult, CCSmult_p, siglev, 0);

					if (neighbors>0) {
						BenjHoch_symmatr (M, CMI,     CMI_p,     siglev, 0);
						BenjHoch_symmatr (M, CMImult, CMImult_p, siglev, 0);
					}
				}

// Holm-Bonferroni / family-wise error rate correction

				if (sw_corr==2) {
					HolmBonf_symmatr (M, CCS,     CCS_p,     siglev, 0);
					HolmBonf_symmatr (M, CCSmult, CCSmult_p, siglev, 0);

					if (neighbors>0) {
						HolmBonf_symmatr (M, CMI,     CMI_p,     siglev, 0);
						HolmBonf_symmatr (M, CMImult, CMImult_p, siglev, 0);
					}
				}

// Bonferroni correction (expected to produce many false negatives!)

				if (sw_corr==3) {
					Bonferroni_symmatr (M, CCS,     CCS_p,     siglev, 0);
					Bonferroni_symmatr (M, CCSmult, CCSmult_p, siglev, 0);

					if (neighbors>0) {
						Bonferroni_symmatr (M, CMI,     CMI_p,     siglev, 0);
						Bonferroni_symmatr (M, CMImult, CMImult_p, siglev, 0);
					}
				}

// *** output of matrices TCS, CCS and RCS ***

				for (i=1;i<=sw_time_out;i++)
					fprintf(fout_matrix_TCS, "  %14.6e    ", t);
				for (i=1;i<=M;i++) {
					for (j=1;j<=M;j++)
						fprintf(fout_matrix_TCS, "  %12.4e", TCS[i][j]);
				}
				fprintf(fout_matrix_TCS, "\n");
				fflush (fout_matrix_TCS);

				for (i=1;i<=sw_time_out;i++)
					fprintf(fout_matrix_RCS, "  %14.6e    ", t);
				for (i=1;i<=M;i++) {
					for (j=1;j<=M;j++)
						fprintf(fout_matrix_RCS, "  %12.4e", RCS[i][j]);
				}
				fprintf(fout_matrix_RCS, "\n");
				fflush (fout_matrix_RCS);

				for (i=1;i<=sw_time_out;i++)
					fprintf(fout_matrix_CCS, "  %14.6e    ", t);
				for (i=1;i<=M;i++) {
					for (j=1;j<=M;j++)
						fprintf(fout_matrix_CCS, "  %12.4e", CCS[i][j]);
				}
				fprintf(fout_matrix_CCS, "\n");
				fflush (fout_matrix_CCS);

				if (neighbors>0) {
					for (i=1;i<=sw_time_out;i++)
						fprintf(fout_matrix_multRCS, "  %14.6e    ", t);
					for (i=1;i<=M;i++) {
						for (j=1;j<=M;j++)
							fprintf(fout_matrix_multRCS, "  %12.4e", RCSmult[i][j]);
					}
					fprintf(fout_matrix_multRCS, "\n");
					fflush (fout_matrix_multRCS);

					for (i=1;i<=sw_time_out;i++)
						fprintf(fout_matrix_multCCS, "  %14.6e    ", t);
					for (i=1;i<=M;i++) {
						for (j=1;j<=M;j++)
							fprintf(fout_matrix_multCCS, "  %12.4e", CCSmult[i][j]);
					}
					fprintf(fout_matrix_multCCS, "\n");
					fflush (fout_matrix_multCCS);

// *** output of matrices TMI, CMI and RMI ***

					for (i=1;i<=sw_time_out;i++)
						fprintf(fout_matrix_TMI, "  %14.6e    ", t);
					for (i=1;i<=M;i++) {
						for (j=1;j<=M;j++)
							fprintf(fout_matrix_TMI, "  %12.4e", TMI[i][j]);
					}
					fprintf(fout_matrix_TMI, "\n");
					fflush (fout_matrix_TMI);

					for (i=1;i<=sw_time_out;i++)
						fprintf(fout_matrix_RMI, "  %14.6e    ", t);
					for (i=1;i<=M;i++) {
						for (j=1;j<=M;j++)
							fprintf(fout_matrix_RMI, "  %12.4e", RMI[i][j]);
					}
					fprintf(fout_matrix_RMI, "\n");
					fflush (fout_matrix_RMI);

					for (i=1;i<=sw_time_out;i++)
						fprintf(fout_matrix_CMI, "  %14.6e    ", t);
					for (i=1;i<=M;i++) {
						for (j=1;j<=M;j++)
							fprintf(fout_matrix_CMI, "  %12.4e", CMI[i][j]);
					}
					fprintf(fout_matrix_CMI, "\n");
					fflush (fout_matrix_CMI);

					for (i=1;i<=sw_time_out;i++)
						fprintf(fout_matrix_multRMI, "  %14.6e    ", t);
					for (i=1;i<=M;i++) {
						for (j=1;j<=M;j++)
							fprintf(fout_matrix_multRMI, "  %12.4e", RMImult[i][j]);
					}
					fprintf(fout_matrix_multRMI, "\n");
					fflush (fout_matrix_multRMI);

					for (i=1;i<=sw_time_out;i++)
						fprintf(fout_matrix_multCMI, "  %14.6e    ", t);
					for (i=1;i<=M;i++) {
						for (j=1;j<=M;j++)
							fprintf(fout_matrix_multCMI, "  %12.4e", CMImult[i][j]);
					}
					fprintf(fout_matrix_multCMI, "\n");
					fflush (fout_matrix_multCMI);
				}
			}
		}
	}

// close output files

	fclose(fin_series);

	fclose(fout_matrix_TCS);
	fclose(fout_matrix_CCS);
	fclose(fout_matrix_RCS);

	if (neighbors>0) {
		fclose(fout_matrix_multCCS);
		fclose(fout_matrix_multRCS);

		fclose(fout_matrix_TMI);
		fclose(fout_matrix_CMI);
		fclose(fout_matrix_RMI);
		fclose(fout_matrix_multCMI);
		fclose(fout_matrix_multRMI);
	}

// deallocate arrays

	free_ivector         (sel,1,columns);
	free_cvector        (line,0,columns*30);

	free_vector                (Xl,1,M);
	free_vector            (Xl_old,1,M);

	free_matrix (X_read,0,T_surr-1,1,M);
	free_matrix (X_surr,0,T_surr-1,1,M);

	free_matrix           (X,0,T-1,1,M);
	free_matrix      (Xtilde,0,T-1,1,M);

	free_f3tensor (X_IFS, 	  0,T_surr-1,1,M,1,N_surr);
	free_matrix   (X_IFSsing, 0,T_surr-1,1,M);
	free_matrix   (C_IFS,                1,M,1,M);
	free_matrix   (MI_IFS,               1,M,1,M);
	free_imatrix  (iter_uni,             1,M,1,N_surr);
	free_matrix   (accur_uni,            1,M,1,N_surr);

	free_f3tensor (X_MFS, 	  0,T_surr-1,1,M,1,N_surr);
	free_matrix   (X_MFSsing, 0,T_surr-1,1,M);
	free_matrix   (C_MFS,                1,M,1,M);
	free_matrix   (MI_MFS,               1,M,1,M);
	free_ivector  (iter_mult,                1,N_surr);
	free_matrix   (accur_mult,           1,M,1,N_surr);

	free_matrix    (C,         1,M,1,M);
	free_f3tensor  (all_C,     1,M,1,M,1,N_ens);
	free_f3tensor  (all_Cabs,  1,M,1,M,1,N_ens);
	free_vector    (sort_C,            1,N_ens);
	free_vector    (sort_Cabs,         1,N_ens);
	free_matrix    (C_med,     1,M,1,M);
	free_matrix    (Cabs_med,  1,M,1,M);

	free_vector    (vec_med,           1,N_ens);

	free_matrix    (C_surr,    1,M,1,M);
	free_f3tensor  (all_C_uni,1,M,1,M,1,N_ens*N_surr);
	free_f3tensor  (all_Cabs_uni,1,M,1,M,1,N_ens*N_surr);
	free_f3tensor  (all_C_mult,1,M,1,M,1,N_ens*N_surr);
	free_f3tensor  (all_Cabs_mult,1,M,1,M,1,N_ens*N_surr);
	free_vector    (sort_C_surr,       1,N_ens*N_surr);
	free_vector    (sort_Cabs_uni,     1,N_ens*N_surr);
	free_vector    (sort_Cabs_mult,    1,N_ens*N_surr);
	free_matrix    (C_uni_med,1,M,1,M);
	free_matrix    (Cabs_uni_med,1,M,1,M);
	free_matrix    (C_mult_med,1,M,1,M);
	free_matrix    (Cabs_mult_med,1,M,1,M);

	free_matrix    (TCS,        1,M,1,M);
	free_matrix    (CCS,        1,M,1,M);
	free_matrix    (RCS,        1,M,1,M);
	free_matrix    (CCSmult,    1,M,1,M);
	free_matrix    (RCSmult,    1,M,1,M);

	free_matrix    (CCS_p,      1,M,1,M);
	free_matrix    (CCSmult_p,  1,M,1,M);

	if (neighbors>0) {
		free_matrix    (MI,         1,M,1,M);
		free_f3tensor  (all_MI,     1,M,1,M,1,N_ens);
		free_f3tensor  (all_MIabs,  1,M,1,M,1,N_ens);
		free_vector    (sort_MI,            1,N_ens);
		free_vector    (sort_MIabs,         1,N_ens);
		free_matrix    (MI_med,     1,M,1,M);
		free_matrix    (MIabs_med,  1,M,1,M);

		free_matrix    (MI_surr,    1,M,1,M);
		free_f3tensor  (all_MI_uni,1,M,1,M,1,N_ens*N_surr);
		free_f3tensor  (all_MIabs_uni,1,M,1,M,1,N_ens*N_surr);
		free_f3tensor  (all_MI_mult,1,M,1,M,1,N_ens*N_surr);
		free_f3tensor  (all_MIabs_mult,1,M,1,M,1,N_ens*N_surr);
		free_vector    (sort_MI_surr,       1,N_ens*N_surr);
		free_vector    (sort_MIabs_uni,     1,N_ens*N_surr);
		free_vector    (sort_MIabs_mult,    1,N_ens*N_surr);
		free_matrix    (MI_uni_med,1,M,1,M);
		free_matrix    (MIabs_uni_med,1,M,1,M);
		free_matrix    (MI_mult_med,1,M,1,M);
		free_matrix    (MIabs_mult_med,1,M,1,M);

		free_matrix    (TMI,         1,M,1,M);
		free_matrix    (CMI,         1,M,1,M);
		free_matrix    (RMI,         1,M,1,M);
		free_matrix    (CMImult,     1,M,1,M);
		free_matrix    (RMImult,     1,M,1,M);

		free_matrix    (CMI_p,       1,M,1,M);
		free_matrix    (CMImult_p,   1,M,1,M);
	}

    return 0;
}

