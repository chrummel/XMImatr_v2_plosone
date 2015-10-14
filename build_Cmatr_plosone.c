#include <stdio.h>
#include <stdlib.h>
#include <math.h>



/************************************************************************************/

/*
 * function builds equal-time (zero-lag) covariance and cross-correlation matrices
 * if time series is normalized cross-correlation matrix and
 * otherwise covariance matrix results
 *
 * power mapping with exponent q is possible
 * reference:
 * P. Andersson, A. Oberg, T. Guhr
 * Power mapping and noise reduction for financial correlations
 * Acta Physica Polonica B36, 2611-2619 (2005)
 */

/*
 * written by Christian Rummel, 2006/09/24, crummel@web.de
 */

/************************************************************************************/

void build_C (int M, int T, double q, double **X, double **C)
{
// data types

	int i, j, n;

// routine

	for (i=1;i<=M;i++) {
		for (j=i;j<=M;j++) {
			C[i][j] = 0.0;
			for (n=0;n<T;n++)
				C[i][j] = C[i][j] + X[n][i]*X[n][j];
			C[i][j] = C[i][j]/T;

			if (q!=1.0)
				C[i][j] = C[i][j]/fabs(C[i][j]) * pow(fabs(C[i][j]),q);		// power mapping

			C[j][i] = C[i][j];
		}
	}
}

