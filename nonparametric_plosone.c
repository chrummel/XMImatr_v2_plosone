#include <stdio.h>
#include <stdlib.h>
#include <math.h>



/************************************************************************************/

/*
 * function for non-parametric U-test by Mann, Whitney and Wilcoxon
 * for equal medians of two distributions
 * significance of U is estimated from normal distribution
 * if distributions are too small, significance 10.0 is returned
 * see Siegel, Nichtparametrische statistische Methoden
 * rho is Herrenstein's statistic for the overlap
 */

/*
 * written by Christian Rummel, 2008/10/15, crummel@web.de
 */

/************************************************************************************/

void Utest (int N1, double distrib1[], int N2, double distrib2[], int *U, double *rho, double *p_normal)
{
//	N1									// size of    1st distribution
//	N2									// size of    2nd distribution
//	distrib1							// vector for 1st distribution
//	distrib2							// vector for 2nd distribution
//	U									// U statistic
//	rho									// Herrenstein's statistic
//	p_normal							// normal approximation to significance of U

// data types

	int n, m;							// indices

	double *distrib_joint; 				// vector for joint distribution
	double *labels;						// vector for labels
	double *ranks;						// vector of ranks

	int small, large;					// indices for borders of ties
	double ave_rank;					// average rank

	double UU;							// U statistic
	double R1, R2;						// rank sum of 1st and 2nd distribution
	double mu_U, sigma_U, z;			// for calculation of normal distributed variant

// prototypes for numerical recipes

	double    *vector            (long, long);
	void  free_vector (double *,  long, long);

	void   sort2 (unsigned long, double[], double[]);
	double erf   (double);

// allocate vectors

	distrib_joint = vector (1,N1+N2);
	labels        = vector (1,N1+N2);
	ranks         = vector (1,N1+N2);

// fill vectors for joint distribution and labels

	for (n=1;n<=N1;n++) {
		distrib_joint[n] = distrib1[n];
		labels       [n] = 1.0;
	}
	for (n=1;n<=N2;n++) {
		distrib_joint[N1+n] = distrib2[n];
		labels       [N1+n] = 2.0;
	}

// sort joint distribution

	sort2 (N1+N2, distrib_joint, labels);

// calculate ranks assigning average rank to ties

	n     = 1;
	small = 1;
	while (n<N1+N2) {
		if (distrib_joint[n+1]>distrib_joint[n]) {
			large = n;

			ave_rank = 0.0;
			for (m=small;m<=large;m++)
				ave_rank = ave_rank + 1.0*m;
			ave_rank = ave_rank / (large-small+1);

			for (m=small;m<=large;m++)
				ranks[m] = ave_rank;

			small = large + 1;
		}

		n++;
	}

	large    = N1 + N2;
	ave_rank = 0.0;
	for (m=small;m<=large;m++)
		ave_rank = ave_rank + 1.0*m;
	ave_rank = ave_rank / (large-small+1);

	for (m=small;m<=large;m++)
		ranks[m] = ave_rank;

// calculate rank sum R1 of 1st distribution and statistic U

	R1 = 0.0;
	R2 = 0.0;
	for (n=1;n<=N1+N2;n++) {
		if ((N1<=N2) && (labels[n]==1.0))
			R1 = R1 + ranks[n];
		if ((N1>N2)  && (labels[n]==2.0))
			R2 = R2 + ranks[n];
	}

	if (N1<=N2)
		UU = N1*N2 + N1*(N1+1)/2 - R1;
	else
		UU = N1*N2 + N2*(N2+1)/2 - R2;

	*U = (int)(UU);

// compute Herrenstein's statistic for the overlap

	*rho = UU/(N1*N2);

// for large enough distributions calculate significance from normal distribution

	if ((N1>20) || (N2>20)) {
		mu_U    =       0.5*N1*N2;
		sigma_U = sqrt (1.0*N1*N2*(N1+N2+1)/12.0);

		z = (UU - mu_U) / sigma_U;

		if (z>=0.0)
			*p_normal = 0.5*(1.0-erf(sqrt(0.5)*z));
		else
			*p_normal = 0.5*(1.0+erf(sqrt(0.5)*z));
	}
	else
		*p_normal = 10.0;

// deallocate vectors

	free_vector (distrib_joint, 1,N1+N2);
	free_vector (labels,        1,N1+N2);
	free_vector (ranks,         1,N1+N2);
}


