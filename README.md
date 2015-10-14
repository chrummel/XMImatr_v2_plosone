# XMImatr_v2
C codes for calculation of surrogate corrected linear and nonlinear interrelation matrices



WARRANTY AND COPYRIGHT:

This toolbox of C programs is made publicly available WITHOUT ANY WARRANTY; without even the implied warranty of merchantability or warranty of fitness for a particular purpose. If you decide to use the toolbox, you do so AT YOUR OWN RISK!

The copyright to each of the distributed programs is held by 

Christian Rummel 
Support Center for Advanced Neuroimaging (SCAN)
University Institute for Diagnostic and Interventional Neuroradiology
Inselspital, 3010 Bern, Switzerland
crummel@web.de



PURPOSE:

The toolbox has been used to generate linear and nonlinear interrelation matrices of intracranial EEG time series for the following publications:

C. Rummel, M. Müller, G. Baier, F. Amor, K. Schindler
Analyzing spatio-temporal patterns of genuine cross-correlations
J. Neurosci. Meth. 191, 94-100 (2010)

C. Rummel, E. Abela, M. Müller, M. Hauf, O. Scheidegger, R. Wiest, K. Schindler
Uniform approach to linear and nonlinear interrelation patterns in multivariate time series
Phys. Rev. E83, 066215 (2011)

C. Rummel, E. Abela, R.G. Andrzejak, M. Hauf, C. Pollo, M. Müller, C. Weisstanner, R. Wiest, K. Schindler
Resected brain tissue, seizure onset zone and quantitative EEG measures:
Towards prediction of post-surgical seizure control
PLoS ONE accepted (2015)

PLEASE CITE these publications if you use these programs! 



CONTENT:

README.txt
build_Cmatr_plosone.c
build_MImatr2_plosone.cpp
multcomp_plosone.c
nonparametric_plosone.c
surrogates_plosone.c
XMImatr_v2_plosone.c
inputs/series_raw_plosone.inp
inputs/XMImatr_v2_plosone.inp



DEPENDENCIES:

The toolbox uses the following third party software packages, which are NOT part of this distribution:

Numerical Recipes in C: The Art of Scientific Computing
by W.H. Press, B.P. Flannery, S.A. Teukolsky, W.T. Vetterling
Cambridge University Press, 2 edition (1992)

The routine mir_xnyn.C is part of the MILCA toolbox by
by S. Astakhov, P. Grassberger, A. Kraskov, H. Stögbauer. 
It can be downloaded at https://www.ucl.ac.uk/ion/departments/sobell/Research/RLemon/MILCA/MILCA (status 2015/10/12). 



USAGE:

Satisfy these dependencies and compile the source code as C++ project. The math library must be linked (-lm). 

The following directories must be provided:
- ../inputs: This directory must contain the files XMImatr_v2_plosone.inp (parameters for analysis) and series_raw_plosone.inp (parameters for input time series). 
- ../data: This directory must contain the input time series series_raw.dat (ASCII file, M data channels as columns). All output matrices is written to this directory (all MxM matrix elements in a single column). 


