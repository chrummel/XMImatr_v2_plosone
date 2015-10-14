#include <stdio.h>
#include <math.h>

#define pi				3.141592653589793



/********************************************************************************/

/*
 * SUBROUTINE
 * for generation of a set of Independent Fourier Based Surrogates (IFS)
 * for a multivariate time series using the
 * Iterated Amplitude Adjusted Fourier Surrogates technique (IAAFT)
 *
 * minimization of end-to-end mismatch is not implemented
 * signals should be NORMALIZED to zero mean and unit variance
 *
 * references:
 *
 * T. Schreiber and A. Schmitz
 * Improved Surrogate Data for Nonlinearity Tests
 * Phys. Rev. Lett. 77, 635 (1996)
 *
 * T. Schreiber and A. Schmitz
 * Surrogate Time Series
 * Physica D142, 346 (2000)
 */

/*
 * written by Christian Rummel, 2007/10/30, crummel@web.de
 */

/********************************************************************************/

void IFS2 (int M, int T, double f_samp, int N_surr, double **X, double ***X_IFS, long seed, int spectrum, int maxiter, double Delta_accur, int **iterat, double **accur)
{
/* arguments:
 * int M								channel number
 * int T								length of time series, must be 2^n
 * double f_samp						sampling frequency in Hz
 * N_surr								number of surrogates to be generated
 * double **X							data matrix for time series, M channels, T samples
 * double ***X_IFS						data cube for IAAFT surrogates
 * long seed							seed for random number generator
 * int spectrum							switch for stopping with identical amplitude distribution (==0) or power spectrum (==1)
 * int maxiter							maximal number of iterations
 * double Delta_accur					stopping criterium: drop of channel-wise discrepancy in amplitudes
 * int **iterat							matrix for number of needed iterations, M channels, N_surr surrogates
 * double **accur						matrix for final accuracies, M channels, N_surr surrogates
 */

// data types

	int i, t, n;						// indices
	int t_unsort, t_0, t_ex;			// time indices
	int iter;							// counter for iterations

	double discrep_iter;				// relative discrepancy accepted on actual iteration step

	double *X_sort, **X_rank;			// rank ordered copy of multivariate time series
	double *ind_unsort;					// vectors for original indices before and after sorting
	double **X_shuffled;				// matrix for shuffled time series
	double **Y;							// matrix for iterated multivariate time series with exact amplitude distributions, complex!
	double *Y_old;						// amplitudes on previous iteration step
	double *Y_uni, *Y_phaserand;		// vectors for univariate time series and its phase randomized copy
	double **Z;							// matrix for iterated multivariate time series with exact Fourier spectra, complex!
	double Delta, Delta1, norm;			// accuracy measure between original data and surrogate

	double **spect, **phase;			// matrices for spectral amplitudes and phases of multivariate time series

// prototypes for numerical recipes routines

	double    *vector (long, long);
	double   **matrix (long, long, long, long);
	void  free_vector (double *, long, long);
	void  free_matrix (double **, long, long, long, long);

	void   sort   (unsigned long, double[]);
	void   sort2  (unsigned long, double[], double[]);
	double ran1   (long *);
	double gasdev (long *);

	void    four1 (double[], unsigned long, int);

// other prototypes

	void Fourier (int, int, int, double **, double **, double **);
	void shuffle (int, int,      double **, double **, long);

// allocate arrays

	spect       = matrix (1,T,1,M);
	phase       = matrix (1,T,1,M);

	Y_uni       = vector (1,2*T);
	Y_phaserand = vector (1,2*T);

	X_sort      = vector (1,T);
	Y_old       = vector (1,T);
	ind_unsort  = vector (1,T);
	X_rank      = matrix (0,T-1,1,M);
	X_shuffled  = matrix (0,T-1,1,M);

// calculate Fourier spectrum of original multivariate time series

	Fourier (0,M,T,X,spect,phase);

// produce rank ordered copy of original multivariate time series for update step 2

	for (i=1;i<=M;i++) {
		for (t=0;t<T;t++)
			X_sort[t+1] = X[t][i];

		sort (T, X_sort);

		for (t=0;t<T;t++)
			X_rank[t][i] = X_sort[t+1];
	}

// loop for generation of surrogate set

	for (n=1;n<=N_surr;n++) {

// start with randomly shuffled multivariate time series
// all channels are shuffled independently

		shuffle (M,T, X,X_shuffled, seed-n);

		Y = matrix (1,2*T,1,M);
		Z = matrix (1,2*T,1,M);

		for (i=1;i<=M;i++) {
			for (t=1;t<=T;t++) {
				Y[2*t-1][i] = X_shuffled[t-1][i];	// real part goes to odd samples
				Y[2*t]  [i] = 0.0;					// zero imaginary part goes to even samples
			}
		}

// loop over independent channels

		for (i=1;i<=M;i++) {
			for (t=1;t<=T;t++)
				Y_old[t] = X_shuffled[t-1][i];

// iterate until convergence or maximal number of iterations is reached
// start with nominally accepted relative discrepancy and increase it by 10% if convergence cannot be reached

			iter = 0;
			while (iter<maxiter) {
				iter++;

// update step 1: ensure correct Fourier spectrum of data matrix Z

				for (t=1;t<=2*T;t++)
					Y_uni[t] = Y[t][i];

				four1 (Y_uni,T,1);

				for (t=1;t<=T;t++) {
					Y_phaserand[2*t-1] = spect[t][i] * cos(atan2(Y_uni[2*t],Y_uni[2*t-1]));
					Y_phaserand[2*t]   = spect[t][i] * sin(atan2(Y_uni[2*t],Y_uni[2*t-1]));
				}

				four1 (Y_phaserand,T,-1);

				for (t=1;t<=2*T;t++)
					Z[t][i] = Y_phaserand[t];

// update step 2: rescale to correct distribution of values of data matrix Y (real part, imaginary part should be zero)

				for (t=1;t<=T;t++) {
					X_sort    [t] = Z[2*t-1][i];
					ind_unsort[t] = 1.0*t;
				}

				sort2 (T, X_sort, ind_unsort);

				for (t=1;t<=T;t++) {
					t_unsort = (int)(ind_unsort[t]);

					Y[2*t_unsort-1][i] = X_rank[t-1][i];
					Y[2*t_unsort]  [i] = 0.0;
				}

// check convergence in amplitudes of real part

				Delta = 0.0;
				for (t=1;t<=T;t++)
					Delta = Delta + (Y_old[t]-Y[2*t-1][i])*(Y_old[t]-Y[2*t-1][i]);
				Delta = sqrt(Delta);

				if (iter==1)
					Delta1 = Delta;
				if (Delta/Delta1<Delta_accur)
					break;

				for (t=1;t<=T;t++)
					Y_old[t] = Y[2*t-1][i];
			}

			if (iter==maxiter)
				printf("\n    IFS2: iteration suspended after %6i steps for channel %4i, accuracy may be low!\n", maxiter, i);

// calculate finally reached accuracies
// if exact amplitude distribution is desired, accuracy is calculated from spectral deviations
// if exact power spectrum is desired, accuracy is calculated from amplitude deviations

			if (spectrum==0) {
				Delta = 0.0;
				norm  = 0.0;

				for (t=1;t<=2*T;t++)
					Y_uni[t] = Y[t][i];

				four1 (Y_uni,T,1);

				for (t=1;t<=T;t++) {
					Delta = Delta + (spect[t][i] - sqrt(Y_uni[2*t-1]*Y_uni[2*t-1]+Y_uni[2*t]*Y_uni[2*t])/T)*(spect[t][i] - sqrt(Y_uni[2*t-1]*Y_uni[2*t-1]+Y_uni[2*t]*Y_uni[2*t])/T);
					norm  = norm  +  spect[t][i]                                                           * spect[t][i];
				}
			}

			if (spectrum==1) {
				Delta = 0.0;
				norm  = 0.0;

				for (t=1;t<=T;t++)
					X_sort[t] = Z[2*t-1][i];

				sort (T, X_sort);

				for (t=1;t<=T;t++) {
					Delta = Delta + (X_rank[t-1][i]-X_sort[t])*(X_rank[t-1][i]-X_sort[t]);
					norm  = norm  +  X_rank[t-1][i]           * X_rank[t-1][i];
				}
			}

			iterat[i][n] = iter;
			accur [i][n] = sqrt(Delta/norm);
		}

// write surrogate to output array
// randomly permute time series about t_0 data points keeping the periodicity

		for (i=1;i<=M;i++) {
			t_0 = (int)(T*ran1(&seed));

			for (t=1;t<=T;t++) {
				if (spectrum==0) {
					if (t+t_0<=T)
						X_IFS[t-1][i][n] = Y[2*t+2*t_0-1]    [i];
					else
						X_IFS[t-1][i][n] = Y[2*t+2*t_0-1-2*T][i];
				}
				else if (spectrum==1) {
					if (t+t_0<=T)
						X_IFS[t-1][i][n] = Z[2*t+2*t_0-1]    [i];
					else
						X_IFS[t-1][i][n] = Z[2*t+2*t_0-1-2*T][i];
				}
			}
		}

		free_matrix (Y, 1,2*T,1,M);
		free_matrix (Z, 1,2*T,1,M);
	}

// deallocate arrays

	free_matrix (spect,       1,T,  1,M);
	free_matrix (phase,       1,T,  1,M);

	free_vector (Y_uni,       1,2*T);
	free_vector (Y_phaserand, 1,2*T);

	free_vector (X_sort,      1,T);
	free_vector (Y_old,       1,T);
	free_vector (ind_unsort,  1,T);
	free_matrix (X_rank,      0,T-1,1,M);
	free_matrix (X_shuffled,  0,T-1,1,M);
}



/********************************************************************************/

/*
 * SUBROUTINE
 * for generation of a set of Multivariate Fourier Based Surrogates (MFS)
 * for a multivariate time series using the
 * Iterated Amplitude Adjusted Fourier Surrogates technique (IAAFT)
 *
 * minimization of end-to-end mismatch is not implemented
 * signals should be NORMALIZED to zero mean and unit variance
 * if convergence cannot be reached, iterations are restarted and
 * accepted discrepancies are slightly increased
 *
 * reference:
 * T. Schreiber and A. Schmitz
 * Surrogate Time Series
 * Physica D142, 346 (2000)
 */

/*
 * written by Christian Rummel, 2010/04/16, crummel@web.de
 */

/********************************************************************************/

void MFS2 (int M, int T, double f_samp, int N_surr, double **X, double ***X_MFS, long seed, int spectrum, int maxiter, double Delta_accur, int iterat[], double **accur)
{
/* arguments:
 * int M								channel number
 * int T								length of time series, must be 2^n
 * double f_samp						sampling frequency in Hz
 * N_surr								number of surrogates to be generated
 * double **X							data matrix for time series, M channels, T samples
 * double ***X_MFS						data cube for IAAFT surrogates
 * long seed							seed for random number generator
 * int spectrum							switch for stopping with identical amplitude distribution (==0) or power spectrum (==1)
 * int maxiter							maximal number of iterations
 * double Delta_accur					stopping criterium: drop of maximal discrepancy in amplitudes
 * int *iterat							vector for number of needed iterations, N_surr surrogates
 * double **accur						matrix for final accuracies, M channels, N_surr surrogates
 */

// data types

	int i, t, n;						// indices
	int t_unsort, t_ex;					// time indices
	int iter;							// counter for iterations

	double discrep_iter;				// discrepancy accepted on actual iteration step
	double Delta_crit;					// critical value of discrepancy used in decisions

	double *X_sort, **X_rank;			// rank ordered copy of multivariate time series
	double *ind_unsort;					// vectors for original indices before and after sorting
	double **X_shuffled;				// matrix for shuffled time series
	double **Y;							// matrix for iterated multivariate time series with exact amplitude distributions, complex!
	double **Y_old;						// amplitudes on previous iteration step
	double *Y_uni, *Y_phaserand;		// vectors for univariate time series and its phase randomized copy
	double **Z;							// matrix for iterated multivariate time series with exact Fourier spectra, complex!
	double Delta, Delta1, norm;			// accuracy measure between original data and surrogate

	double **spect, **phase;			// matrices for spectral amplitudes and phases of multivariate time series
	double **psi;						// matrix for phases
	double h1, h2, alph;
	double *num, *den, *alpha;			// vector for optimal phase shifts

// prototypes for numerical recipes routines

	double    *vector (long, long);
	double   **matrix (long, long, long, long);
	void  free_vector (double *, long, long);
	void  free_matrix (double **, long, long, long, long);

	void   sort   (unsigned long, double[]);
	void   sort2  (unsigned long, double[], double[]);
	double ran1   (long *);
	double gasdev (long *);

	void    four1 (double[], unsigned long, int);

// other prototypes

	void Fourier  (int, int, int, double **, double **, double **);
	void shuffle  (int, int,      double **, double **, long);
	void shuffle1 (int, int,      double **, double **, long);

// allocate arrays

	spect       = matrix (1,T,1,M);
	phase       = matrix (1,T,1,M);

	Y_uni       = vector (1,2*T);
	Y_phaserand = vector (1,2*T);

	X_sort      = vector (1,T);
	Y_old       = matrix (1,T,1,M);
	ind_unsort  = vector (1,T);
	X_rank      = matrix (0,T-1,1,M);
	X_shuffled  = matrix (0,T-1,1,M);

	psi         = matrix (1,T,1,M);
	den         = vector (1,T);
	num         = vector (1,T);
	alpha       = vector (1,T);

// calculate Fourier spectrum of original multivariate time series

	Fourier (0,M,T,X,spect,phase);

// produce rank ordered copy of original multivariate time series for update step 2

	for (i=1;i<=M;i++) {
		for (t=0;t<T;t++)
			X_sort[t+1] = X[t][i];

		sort (T, X_sort);

		for (t=0;t<T;t++)
			X_rank[t][i] = X_sort[t+1];
	}

// loop for generation of surrogate set

	for (n=1;n<=N_surr;n++) {

// start with randomly shuffled multivariate time series
// all channels are shuffled simultaneously
// this is THE FIRST of three places where the multivariate version is essentially different from the univariate version

		shuffle1 (M,T, X,X_shuffled, seed-n);

		Y = matrix (1,2*T,1,M);
		Z = matrix (1,2*T,1,M);

		for (i=1;i<=M;i++) {
			for (t=1;t<=T;t++) {
				Y[2*t-1][i] = X_shuffled[t-1][i];	// real part goes to odd samples
				Y[2*t]  [i] = 0.0;					// zero imaginary part goes to even samples

				Y_old[t][i] = X_shuffled[t-1][i];
			}
		}

// iterate until convergence or maximal number of iterations is reached
// start with nominally accepted relative discrepancy and increase it by 10% if convergence cannot be reached

		iter = 0;
		while (iter<maxiter) {
			iter++;

// calculate optimal phase shifts
// this is THE SECOND of three places where the multivariate version is essentially different from the univariate version

			for (t=1;t<=T;t++) {
				num[t] = 0.0;
				den[t] = 0.0;
			}

			for (i=1;i<=M;i++) {
				for (t=1;t<=2*T;t++)
					Y_uni[t] = Y[t][i];

				four1 (Y_uni,T,1);

				for (t=1;t<=T;t++) {
					psi[t][i] = atan2(Y_uni[2*t],Y_uni[2*t-1]);
					num[t]    = num[t] + sin(psi[t][i]-phase[t][i]);
					den[t]    = den[t] + cos(psi[t][i]-phase[t][i]);
				}
			}

			for (t=1;t<=T;t++) {
				alph = atan2(num[t],den[t]);

				h1 = 0.0;
				h2 = 0.0;
				for (i=1;i<=M;i++) {
					h1 = h1 + 2.0 - 2.0*cos(alph   -psi[t][i]+phase[t][i]);
					h2 = h2 + 2.0 - 2.0*cos(alph+pi-psi[t][i]+phase[t][i]);
				}

				if (h1<=h2)
					alpha[t] = alph;
				else
					alpha[t] = alph + pi;
			}

// loop over univariate time series

			Delta_crit = 0.0;
			for (i=1;i<=M;i++) {

// update step 1: ensure correct Fourier spectrum of data matrix Z
// this is THE THIRD of three places where the multivariate version is essentially different from the univariate version

				for (t=1;t<=T;t++) {
					Y_phaserand[2*t-1] = spect[t][i] * cos(phase[t][i]+alpha[t]);
					Y_phaserand[2*t]   = spect[t][i] * sin(phase[t][i]+alpha[t]);
				}

				four1 (Y_phaserand,T,-1);

				for (t=1;t<=2*T;t++)
					Z[t][i] = Y_phaserand[t];

// update step 2: rescale to correct distribution of values of data matrix Y (real part, imaginary part should be zero)

				for (t=1;t<=T;t++) {
					X_sort    [t] = Z[2*t-1][i];
					ind_unsort[t] = 1.0*t;
				}

				sort2 (T, X_sort, ind_unsort);

				for (t=1;t<=T;t++) {
					t_unsort = (int)(ind_unsort[t]);

					Y[2*t_unsort-1][i] = X_rank[t-1][i];
					Y[2*t_unsort]  [i] = 0.0;
				}

// check convergence in amplitudes of real part

				Delta = 0.0;
				for (t=1;t<=T;t++)
					Delta = Delta + (Y_old[t][i]-Y[2*t-1][i])*(Y_old[t][i]-Y[2*t-1][i]);
				Delta = sqrt(Delta);

				if (Delta>Delta_crit)
					Delta_crit = Delta;

				for (t=1;t<=T;t++)
					Y_old[t][i] = Y[2*t-1][i];
			}

			if (iter==1)
				Delta1 = Delta_crit;
			if (Delta_crit/Delta1<Delta_accur)
				break;
		}

		if (iter==maxiter)
			printf("\n    MFS2: iteration suspended after %6i steps, accuracy may be low!\n", maxiter);

// calculate finally reached accuracies
// if exact amplitude distribution is desired, accuracy is calculated from spectral deviations
// if exact power spectrum is desired, accuracy is calculated from amplitude deviations

		for (i=1;i<=M;i++) {
			if (spectrum==0) {
				Delta = 0.0;
				norm  = 0.0;

				for (t=1;t<=2*T;t++)
					Y_uni[t] = Y[t][i];

				four1 (Y_uni,T,1);

				for (t=1;t<=T;t++) {
					Delta = Delta + (spect[t][i] - sqrt(Y_uni[2*t-1]*Y_uni[2*t-1]+Y_uni[2*t]*Y_uni[2*t])/T)*(spect[t][i] - sqrt(Y_uni[2*t-1]*Y_uni[2*t-1]+Y_uni[2*t]*Y_uni[2*t])/T);
					norm  = norm  +  spect[t][i]                                                           * spect[t][i];
				}
			}

			if (spectrum==1) {
				Delta = 0.0;
				norm  = 0.0;

				for (t=1;t<=T;t++)
					X_sort[t] = Z[2*t-1][i];

				sort (T, X_sort);

				for (t=1;t<=T;t++) {
					Delta = Delta + (X_rank[t-1][i]-X_sort[t])*(X_rank[t-1][i]-X_sort[t]);
					norm  = norm  +  X_rank[t-1][i]           * X_rank[t-1][i];
				}
			}

			accur[i][n] = sqrt(Delta/norm);
		}

		iterat[n] = iter;

// write surrogate to output array

		for (i=1;i<=M;i++) {
			for (t=1;t<=T;t++) {
				if (spectrum==0)
					X_MFS[t-1][i][n] = Y[2*t-1][i];
				else if (spectrum==1)
					X_MFS[t-1][i][n] = Z[2*t-1][i];
			}
		}

		free_matrix (Y, 1,2*T,1,M);
		free_matrix (Z, 1,2*T,1,M);
	}

// clean up

	free_matrix (spect,       1,T,  1,M);
	free_matrix (phase,       1,T,  1,M);

	free_matrix (psi,         1,T,  1,M);
	free_vector (num,         1,T);
	free_vector (den,         1,T);
	free_vector (alpha,       1,T);

	free_vector (Y_uni,       1,2*T);
	free_vector (Y_phaserand, 1,2*T);

	free_vector (X_sort,      1,T);
	free_matrix (Y_old,       1,T,  1,M);
	free_vector (ind_unsort,  1,T);
	free_matrix (X_rank,      0,T-1,1,M);
	free_matrix (X_shuffled,  0,T-1,1,M);
}



/********************************************************************************/

/*
 * SUBROUTINE
 * for simple shuffling of the time line of a multivariate time series
 * all channels are shuffled independently
 */

/*
 * written by Christian Rummel, 2007/10/30, crummel@web.de
 */

/********************************************************************************/

void shuffle (int M, int T, double **X, double **Y, long seed)
{
// data types

	int i, t;							// indices
	int t_ex;							// index of randomly chosen time step
	double Y_ex;						// dummy for exchanging entries

// prototypes for numerical recipes routines

	double    *vector (long, long);
	double   **matrix (long, long, long, long);
	void  free_vector (double *, long, long);
	void  free_matrix (double **, long, long, long, long);

	double ran1   (long *);

// fill output matrix with copy of original time series

	for (i=1;i<=M;i++) {
		for (t=0;t<T;t++)
			Y[t][i] = X[t][i];
	}

// for each time step exchange the matrix element with the one of a randomly chosen time step

	for (t=0;t<T;t++) {
		for (i=1;i<=M;i++) {
			t_ex       = (int)(T*ran1(&seed));

			Y_ex       = Y[t]   [i];
			Y[t]   [i] = Y[t_ex][i];
			Y[t_ex][i] = Y_ex;
		}
	}
}



/********************************************************************************/

/*
 * SUBROUTINE
 * for simple shuffling of the time line of a multivariate time series
 * all channels are shuffled simultaneously
 */

/*
 * written by Christian Rummel, 2010/08/18, crummel@web.de
 */

/********************************************************************************/

void shuffle1 (int M, int T, double **X, double **Y, long seed)
{
// data types

	int i, t;							// indices
	int t_ex;							// index of randomly chosen time step
	double Y_ex;						// dummy for exchanging entries

// prototypes for numerical recipes routines

	double    *vector (long, long);
	double   **matrix (long, long, long, long);
	void  free_vector (double *, long, long);
	void  free_matrix (double **, long, long, long, long);

	double ran1   (long *);

// fill output matrix with copy of original time series

	for (i=1;i<=M;i++) {
		for (t=0;t<T;t++)
			Y[t][i] = X[t][i];
	}

// for each time step exchange the matrix element with the one of a randomly chosen time step

	for (t=0;t<T;t++) {
		t_ex = (int)(T*ran1(&seed));

		for (i=1;i<=M;i++) {
			Y_ex       = Y[t]   [i];
			Y[t]   [i] = Y[t_ex][i];
			Y[t_ex][i] = Y_ex;
		}
	}
}



/********************************************************************************/

/*
 * SUBROUTINE
 * for calculation of amplitudes and phases of a multivariate time series
 * frequencies are arranged as in routine four1 from Numerical Recipes:
 * - positive frequencies from small to large for t=1...T/2
 * - negative frequencies from large to small for t=T/2+1...T
 */

/*
 * written by Christian Rummel, 2007/10/30, crummel@web.de
 */

/********************************************************************************/

void Fourier (int W, int M, int T, double **X, double **spect, double **phase)
{
// data types

	int i, t;							// indices

	double *X_uni;						// vector for univariate time series

// prototypes for numerical recipes routines

	double    *vector (long, long);
	void  free_vector (double *, long, long);

	void        four1 (double[], unsigned long, int);

// loop for time series

	X_uni = vector (1,2*T);

	for (i=1;i<=M;i++) {

// copy univariate time series to vector

		for (t=1;t<=T;t++) {
			if (W==0)
				X_uni[2*t-1] = X[t-1][i];													// square window
			else if (W==1)
				X_uni[2*t-1] = X[t-1][i] * (1.0 - 4.0*(t-T/2)*(t-T/2)/T/T);					// Welch window
			X_uni[2*t] = 0.0;
		}

// Fast Fourier Transform of univariate time series

		four1 (X_uni,T,1);

// calculate Fourier amplitude and phase of univariate time series
// real part in odd and imaginary part in even samples

		for (t=1;t<=T;t++) {
			spect[t][i] = sqrt  (X_uni[2*t-1]*X_uni[2*t-1] + X_uni[2*t]*X_uni[2*t]) / T;
			phase[t][i] = atan2 (X_uni[2*t],X_uni[2*t-1]);
		}
	}

	free_vector (X_uni, 1,2*T);
}



#undef	pi
