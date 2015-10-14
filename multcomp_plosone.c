#include <stdio.h>
#include <stdlib.h>

#include <math.h>



/********************************************************************************/

/*
 * subroutines for significance thresholding
 * and multiple comparison correction of symmetric matrices
*/

/*
 * written by Christian Rummel, crummel@web.de
*/

/********************************************************************************/



/*
 * Bonferroni correction (expected to produce many false negatives!)
 */

void Bonferroni_symmatr (int M, double **matr, double **p, double siglev, int hypo)
{
//	hypo==0 rejects null hypothesis, hypo==1 rejects alternative hypothesis

	int i, j;				// index

// Bonferroni correction

	for (i=1;i<=M;i++) {
		for (j=i+1;j<=M;j++)
			p[i][j] = p[i][j]*M*(M-1)/2;
	}

// delete insignificant elements

	if (hypo==0) {
		for (i=1;i<=M;i++) {
			for (j=i+1;j<=M;j++) {
				if (p[i][j]>=siglev)
					matr[i][j] = 0.0;
			}
		}
	}
	else if (hypo==1) {
		for (i=1;i<=M;i++) {
			for (j=i+1;j<=M;j++) {
				if (p[i][j]<siglev)
					matr[i][j] = 0.0;
			}
		}
	}

	for (i=1;i<=M;i++) {
		for (j=i+1;j<=M;j++) {
			matr[j][i] = matr[i][j];
			if (p[i][j]>1.0)
				p[i][j] = 1.0;
			p[j][i] = p[i][j];
		}
	}
}



/*
 * Holm-Bonferroni / family-wise error rate correction
 */

void HolmBonf_symmatr (int M, double **matr, double **p, double siglev, int hypo)
{
//	hypo==0 rejects null hypothesis, hypo==1 rejects alternative hypothesis

	int i, j, k;			// index
	int N_sig;				// number of significant elements

	double *sort_v;			// vector for sorting
	double *sort_p;			// vector for sorting p-values
	double *sort_i;			// vector for sorting indices
	double *sort_j;			// vector for sorting indices

	double     *vector (long, long);
	void   free_vector (double *, long, long);

	void   sort2 (int, double[], double[]);

// allocate vectors

	sort_v = vector (1,M*(M-1)/2);
	sort_p = vector (1,M*(M-1)/2);
	sort_i = vector (1,M*(M-1)/2);
	sort_j = vector (1,M*(M-1)/2);

// sort p-values and indices

	k = 0;
	for (i=1;i<=M;i++) {
		for (j=i+1;j<=M;j++) {
			k++;
			sort_v[k] = p[i][j];
			sort_i[k] = 1.0*i;
		}
	}
	sort2 (M*(M-1)/2, sort_v, sort_i);

	k = 0;
	for (i=1;i<=M;i++) {
		for (j=i+1;j<=M;j++) {
			k++;
			sort_v[k] = p[i][j];
			sort_j[k] = 1.0*j;
		}
	}
	sort2 (M*(M-1)/2, sort_v, sort_j);

// step down procedure:
// find number of significant elements and correct p-values

	for (k=1;k<=M*(M-1)/2;k++) {
		sort_v[k] = sort_v[k] * (0.5*M*(M-1)-k+1);
		p[(int)(sort_i[k])][(int)(sort_j[k])] = sort_v[k];
	}
	N_sig = 0;
	while ((N_sig<M*(M-1)/2) && (sort_v[N_sig+1]<siglev))
		N_sig++;

// delete insignificant elements

	if (hypo==0) {
		for (k=N_sig+1;k<=M*(M-1)/2;k++)
			matr[(int)(sort_i[k])][(int)(sort_j[k])] = 0.0;
	}
	else if (hypo==1) {
		for (k=1;k<=N_sig;k++)
			matr[(int)(sort_i[k])][(int)(sort_j[k])] = 0.0;
	}

	for (i=1;i<=M;i++) {
		for (j=i+1;j<=M;j++) {
			matr[j][i] = matr[i][j];
			if (p[i][j]>1.0)
				p[i][j] = 1.0;
			p[j][i] = p[i][j];
		}
	}

// deallocate vectors

	free_vector (sort_v, 1,M*(M-1)/2);
	free_vector (sort_p, 1,M*(M-1)/2);
	free_vector (sort_i, 1,M*(M-1)/2);
	free_vector (sort_j, 1,M*(M-1)/2);
}



/*
 * Benjamini-Hochberg / false discovery rate correction
 */

void BenjHoch_symmatr (int M, double **matr, double **p, double siglev, int hypo)
{
//	hypo==0 rejects null hypothesis, hypo==1 rejects alternative hypothesis

	int i, j, k;			// index
	int N_sig;				// number of significant elements

	double *sort_v;			// vector for sorting
	double *sort_p;			// vector for sorting p-values
	double *sort_i;			// vector for sorting indices
	double *sort_j;			// vector for sorting indices

	double    *vector            (long, long);
	void  free_vector (double *,  long, long);

	void   sort2 (int, double[], double[]);

// allocate vectors

	sort_v = vector (1,M*(M-1)/2);
	sort_p = vector (1,M*(M-1)/2);
	sort_i = vector (1,M*(M-1)/2);
	sort_j = vector (1,M*(M-1)/2);

// sort p-values and indices

	k = 0;
	for (i=1;i<=M;i++) {
		for (j=i+1;j<=M;j++) {
			k++;
			sort_v[k] = p[i][j];
			sort_i[k] = 1.0*i;
		}
	}
	sort2 (M*(M-1)/2, sort_v, sort_i);

	k = 0;
	for (i=1;i<=M;i++) {
		for (j=i+1;j<=M;j++) {
			k++;
			sort_v[k] = p[i][j];
			sort_j[k] = 1.0*j;
		}
	}
	sort2 (M*(M-1)/2, sort_v, sort_j);

// step up procedure:
// find number of significant elements and correct p-values

	for (k=1;k<=M*(M-1)/2;k++) {
		sort_v[k] = sort_v[k] * 0.5*M*(M-1)/k;
		p[(int)(sort_i[k])][(int)(sort_j[k])] = sort_v[k];
	}
	N_sig = M*(M-1)/2;
	while ((N_sig>0) && (sort_v[N_sig]>=siglev))
		N_sig--;

// delete insignificant elements

	if (hypo==0) {
		for (k=N_sig+1;k<=M*(M-1)/2;k++)
			matr[(int)(sort_i[k])][(int)(sort_j[k])] = 0.0;
	}
	else if (hypo==1) {
		for (k=1;k<=N_sig;k++)
			matr[(int)(sort_i[k])][(int)(sort_j[k])] = 0.0;
	}

	for (i=1;i<=M;i++) {
		for (j=i+1;j<=M;j++) {
			matr[j][i] = matr[i][j];
			if (p[i][j]>1.0)
				p[i][j] = 1.0;
			p[j][i] = p[i][j];
		}
	}

// deallocate vectors

	free_vector (sort_v, 1,M*(M-1)/2);
	free_vector (sort_p, 1,M*(M-1)/2);
	free_vector (sort_i, 1,M*(M-1)/2);
	free_vector (sort_j, 1,M*(M-1)/2);
}



/*
 * significance check without correction for multiple tests
 */

void uncorr_symmatr (int M, double **matr, double **p, double siglev, int hypo)
{
//	hypo==0 rejects null hypothesis, hypo==1 rejects alternative hypothesis

	int i, j;				// index

// delete insignificant elements

	if (hypo==0) {
		for (i=1;i<=M;i++) {
			for (j=i+1;j<=M;j++) {
				if (p[i][j]>=siglev)
					matr[i][j] = 0.0;
			}
		}
	}
	else if (hypo==1) {
		for (i=1;i<=M;i++) {
			for (j=i+1;j<=M;j++) {
				if (p[i][j]<siglev)
					matr[i][j] = 0.0;
			}
		}
	}

	for (i=1;i<=M;i++) {
		for (j=i+1;j<=M;j++) {
			matr[j][i] = matr[i][j];
			p   [j][i] = p   [i][j];
		}
	}
}


