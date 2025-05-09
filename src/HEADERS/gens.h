#ifndef GENS_H
#define GENS_H

#include "config.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>

// definitions of success and failure
#define SUCCESS GSL_SUCCESS
#define FAILURE !SUCCESS

// fit types
typedef enum { 
  UNWEIGHTED , UNCORRELATED , CORRELATED
} corrtype ;

// effective mass types
typedef enum {
  LOG_EFFMASS , LOGFWD_EFFMASS , LOGBWD_EFFMASS , LOG2_EFFMASS ,
  ACOSH_EFFMASS , ASINH_EFFMASS , ATANH_EFFMASS , ACOSH_ITERATIVE_EFFMASS ,
  ASINH_ITERATIVE_EFFMASS , EVALUE_EFFMASS } effmass_type ;

// fit types
typedef enum {
  ALPHA_D0 , ALPHA_D0_MULTI , ADLERALPHA_D0 , ADLERALPHA_D0_MULTI , HALEXP , EXP , EXP_XINV , COSH , COSH_ASYMM, COSH_PLUSC , EXP_PLUSC , HLBL_CONT , NRQCD_EXP , NRQCD_EXP2 , NOFIT , PADE , PEXP , POLY , PP_AA , PP_AA_WW , PP_AA_WW_R2 , PP_AA_EXP , PPAA , QCORR_BESSEL , QSUSC_SU2 , SINH , TANH , POLES , QSLAB , QSLAB_FIXED , CORNELL , CORNELL_V2 , FVOL1 , FVOL2 , FVOL3 , FVOL4, FVOL5, FVOL6, UDCB_HEAVY , C4C7 , SOL , SOL2, SU2_SHITFIT , SUN_CONT , ZV_EXP , FVOLCC, LARGENB, FVOL_DELTA , TEST
} fittype ;

// time folding types
typedef enum {
  PLUS_PLUS , PLUS_MINUS , MINUS_PLUS , MINUS_MINUS , NOFOLD , NOFOLD_MINUS , TDER , NOFOLD_SWAPT , NOFOLD_MINUS_SWAPT , ZV_SUB 
} fold_type ;

// special channel types
enum { Vi = 123 , Ai = 678 , Tij = 10203 , Tit = 31323 } ;

// what type of data do we use
typedef enum { Raw , JackKnife , BootStrap } resample_type ;

// file type we expect to read
typedef enum { Corr_File , Distribution_File , Fake_File , Flat_File , GLU_Tcorr_File , GLU_File , GLU_Qmoment_File , Adler_File } file_type ;

typedef enum { Adler , Alphas , Beta_crit , Binding_Corr , Correlator , Exceptional , Fit , Fpi_CLS , General , HLBL , HVP , KKops , KK_BK , Nrqcd , PCAC, Pof , Qcorr , Qsusc , Qslab , QslabFix , Ren_Rats , SpinOrbit, StaticPotential , TetraGEVP , TetraGEVP_Fixed , Wflow , Sol , ZV } analysis_type ;

// x-data descriptor
struct x_desc {
  double X ;
  size_t LT ;
  size_t N ;
  size_t M ;
} ;

// little histogram struct
struct histogram {
  double binmid ;
  double width ;
  size_t count ;
} ;

// map structure
struct pmap {
  size_t *p ;  // is the parameter map
  size_t bnd ; // is the simultaneous index
} ;

// priors struct
struct prior {
  double Val ;
  double Err ;
  bool Initialised ;
} ;

// fit function stuff
struct ffunction {
  double *f ; // ( f_i(a) - y_i )
  double **df ; // first derivatives
  double **d2f ; // second derivatives
  double *fparams ; // fit parameters
  double **U ; // linear fit matrix
  const struct prior *Prior ;  
  size_t N ; // length of the fit
  size_t NPARAMS ; // number of fit parameters
  corrtype CORRFIT ; // type of fit
  double chisq ; // chisq 
} ;

// struct containing our statistics
struct resampled {
  double *resampled ;
  double avg ;
  double err_hi ;
  double err_lo ;
  double err ;
  size_t NSAMPLES ;
  resample_type restype ;
} ;

// struct describing our correlation matrix
struct correlation {
  double **W ;
  double Eigenvalue_Tol ;
  bool Divided_Covariance ;
  bool Column_Balanced ;
} ;

// struct containing the data information
struct data_info {
  struct resampled *x ;
  struct resampled *y ;
  struct correlation Cov ;
  size_t *LT ;
  size_t Nboots ;
  size_t *Ndata ;
  size_t Nsim ;
  size_t Ntot ;
  resample_type Restype ;
} ;

// struct for keeping the fit information
struct fit_info {
  corrtype Corrfit ;
  fittype Fitdef ;
  double *Guess ;
  bool Guesses_Initialised ;
  size_t M ;
  struct pmap *map ;
  int (*Minimize) ( void *fdesc ,
		    const void *data ,
		    const double **W ,
		    const double TOL ) ;
  size_t N ;
  size_t Nlogic ;
  size_t Nparam ;
  size_t Nprior ;
  struct prior *Prior ;
  bool *Sims ;
  double Tol ;
} ;

// fit descriptor struct
struct fit_descriptor {
  struct ffunction f ;
  double (*func)( const struct x_desc X , const double *fparams , const size_t Npars ) ;
  void (*F) ( double *f , const void *data , const double *fparams ) ;
  void (*dF) ( double **df , const void *data , const double *fparams ) ;
  void (*d2F) ( double **d2f , const void *data , const double *fparams ) ;
  void (*guesses) ( double *fparams , const struct data_info Data , const struct fit_info Fit ) ;
  void (*linmat) ( double **U , const void *data , const size_t N , const size_t M , const size_t Nlogic ) ;
  const struct prior *Prior ;
  size_t Nparam ; // Number of parameters
  size_t Nlogic ;  // logical Nparameters
  size_t N ;
  size_t M ;
} ;

// data structure in the fits
struct data {
  size_t n ;
  double *x ;
  double *y ;
  size_t *LT ;
  size_t Npars ;
  struct pmap *map ;
  size_t N ;
  size_t M ;
} ;

// uninitialised flag
#define UNINIT_FLAG (123456789)

// tokenize the input file
struct flat_file {
  char *Token ;
  size_t Token_Length ;
  char *Value ;
  size_t Value_Length ;
} ;

// trajectory information
struct traj {
  size_t Begin ;
  size_t Bin ;
  size_t *Dimensions ;
  size_t End ;
  char *FileX ;
  char *FileY ;
  char *RW ;
  double Fit_High ;
  double Fit_Low ;
  size_t Gs ;
  size_t Gk ;
  fold_type Fold ;
  size_t Increment ;
  size_t Nd ;
  double *mom ;
} ;

// simple graph stuff
struct graph {
  char *Name ;
  char *Xaxis ;
  char *Yaxis ;
  size_t Granularity ;
} ;

// input parameters
struct input_params {
  struct data_info Data ;
  analysis_type Analysis ;
  file_type FileType ;
  struct fit_info Fit ;
  struct graph Graph ;
  struct traj *Traj ;
} ;

#endif
