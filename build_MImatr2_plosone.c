/************************************************************************************/

/*
 *  function builds mutual information matrix on basis of MILCA algorithm by
 *  A. Kraskov, H. Stoegbauer, P. Grassberger
 *  Estimating Mutual Information
 *  Phys. Rev. E69, 066138 (2004)
 */

/*
 *  The function mir_xnyn is part of the file miutils.C in the MILCA toolbox
 *  which as third party software is *not* part of this distribution
 *  it can be downloaded at (status 2015/10/12)
 *  https://www.ucl.ac.uk/ion/departments/sobell/Research/RLemon/MILCA/MILCA
 */

/*
 * written by Christian Rummel, 2007/08/16, crummel@web.de
 */

/************************************************************************************/

extern "C" void build_MI (int, int, int, double, double **, double **, long);		// make subroutine callable from C

// prototype for C subroutines

extern "C" {
	double   **matrix (long, long, long, long);
	void  free_matrix (double **, long, long, long, long);

	double gasdev (long *);
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#define eps		1.0e-12



void build_MI (int M, int T, int k, double q, double **Xtilde, double **NMI, long seed)
{
// data types

	int d;
	int i, j, t;

	double **MI;

	double **x;
	double *scal;
	double *min;
	double *max;
	double *psi;

	double mir;
	int BOX1;

// MILCA prototypes

	void mir_xnyn (double **, int, int, int, int, double *, double *, double *);

// fill mutual information matrix without normalization (evaluation of MI adapted from MILCA package)

	MI = matrix (1,M,1,M);

	x=(double**)calloc(2,sizeof(double*));
	for (d=0;d<2;d++) x[d]=(double*)calloc(T,sizeof(double));
	scal=(double*)calloc(2,sizeof(double));
	min=(double*)calloc(2,sizeof(double));
	max=(double*)calloc(2,sizeof(double));
	psi=(double*)calloc(T+1,sizeof(double));

	psi[1]=-(double).57721566490153;
	for (t=1;t<T;t++) psi[t+1]=psi[t]+1/(double)t;
	BOX1=T-5;

	for (i=1;i<=M;i++) {
		for (d=0;d<2;d++) {min[d]=DBL_MAX/2;max[d]=-DBL_MAX/2;}
	    for (t=0;t<T;t++) {
			x[0][t] = Xtilde[t][i] + eps * gasdev(&seed);
			x[1][t] = Xtilde[t][i] + eps * gasdev(&seed);
		}
		for (d=0;d<2;d++) {
			for (t=0;t<T;t++) {
				if (x[d][t]<min[d]) min[d]=x[d][t];
				if (x[d][t]>max[d]) max[d]=x[d][t];
			}
			for (t=0;t<T;t++)
				x[d][t]=x[d][t]-min[d];
		}

		for (d=0;d<2;d++) scal[d]=BOX1/(max[d]-min[d]);
		mir_xnyn (x,1,1,T,k,psi,scal,&mir);
		MI[i][i] = mir;
	}

	for (i=1;i<=M;i++) {
		for (j=i+1;j<=M;j++) {
			for (d=0;d<2;d++) {min[d]=DBL_MAX/2;max[d]=-DBL_MAX/2;}
		    for (t=0;t<T;t++) {
				x[0][t] = Xtilde[t][i] + eps * gasdev(&seed);
				x[1][t] = Xtilde[t][j] + eps * gasdev(&seed);
			}
			for (d=0;d<2;d++) {
				for (t=0;t<T;t++) {
					if (x[d][t]<min[d]) min[d]=x[d][t];
					if (x[d][t]>max[d]) max[d]=x[d][t];
				}
				for (t=0;t<T;t++)
					x[d][t]=x[d][t]-min[d];
			}

			for (d=0;d<2;d++) scal[d]=BOX1/(max[d]-min[d]);
			mir_xnyn (x,1,1,T,k,psi,scal,&mir);
			MI[i][j] = mir;
		}
	}

// output of normalized mutual information matrix using Joe's formula

	for (i=1;i<=M;i++) {
		for (j=i;j<=M;j++) {
			if (MI[i][j]>=0.0)
				NMI[i][j] = sqrt(1.0 - exp(-2.0*MI[i][j]));
			else
				NMI[i][j] = 0.0;
			NMI[j][i] = NMI[i][j];
		}
	}

	free_matrix (MI, 1,M,1,M);

	for (d=0;d<2;d++) free(x[d]); free(x);
	free(scal);
	free(min);free(max);
	free(psi);
}



#undef eps
