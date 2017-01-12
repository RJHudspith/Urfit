#ifndef GENS_H
#define GENS_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

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
  LOG_EFFMASS , LOG2_EFFMASS , ACOSH_EFFMASS , 
  ASINH_EFFMASS , ATANH_EFFMASS , EVALUE_EFFMASS } effmass_type ;

// fit types
typedef enum {
  EXP , COSH , SINH , EXP_PLUSC , PADE , POLY
} fittype ;

// time folding types
typedef enum {
  PLUS_PLUS , PLUS_MINUS , MINUS_PLUS , MINUS_MINUS , NOFOLD 
} foldtype ;

// error types
enum { ERR , HI , LO , AVE } errtype ;

// what type of data do we use
typedef enum { RAWDATA , JACKDATA , BOOTDATA } resample_type ;

// x-data descriptor
struct x_desc {
  double X ;
  size_t LT ;
} ;

// map structure
struct pmap {
  size_t *p ;
} ;

// fit function stuff
struct ffunction {
  double *f ; // ( f_i(a) - y_i )
  double **df ; // first derivatives
  double **d2f ; // second derivatives
  double *fparams ; // fit parameters
  double *prior ; // prior parameters
  double *err_prior ; // prior errors
  size_t N ; // length of the fit
  size_t NPARAMS ; // number of fit parameters
  corrtype CORRFIT ; // type of fit
  double chisq ; // chisq 
} ;

struct fit_descriptor {
  struct ffunction f ;
  double (*func)( const struct x_desc X , const double *fparams , const size_t Npars ) ;
  void (*F) ( double *f , const void *data , const double *fparams ) ;
  void (*dF) ( double **df , const void *data , const double *fparams ) ;
  void (*d2F) ( double **d2f , const void *data , const double *fparams ) ;
  void (*guesses) ( double *fparams , const size_t Nlogic ) ;
  void (*set_priors) ( double *priors , double *err_priors ) ;
  size_t Nparam ; // Number of parameters
  size_t Nlogic ;  // logical Nparameters 
} ;

// priors struct
struct prior {
  double Val ;
  double Err ;
  bool Initialised ;
} ;

// struct for keeping the fit information
struct fit_info {
  corrtype Corrfit ;
  fittype Fitdef ;
  size_t M ;
  struct pmap *map ;
  int (*Minimize) ( struct fit_descriptor *fdesc ,
		    const void *data ,
		    const double **W ,
		    const double TOL ) ;
  size_t N ;
  size_t Nlogic ;
  size_t Nparam ;
  struct prior *Prior ;
  bool *Sims ;
  double Tol ;
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
  size_t LT ;
  size_t Nboots ;
  size_t *Ndata ;
  size_t Nsim ;
  size_t Ntot ;
} ;

// data structure in the fits
struct data {
  size_t n ;
  double *x ;
  double *y ;
  size_t LT ;
  size_t Npars ;
  struct pmap *map ;
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
  char *Filename ;
  double Fit_High ;
  double Fit_Low ;
  size_t Increment ;
  size_t Nd ;
} ;

struct fit_inputs {
  fittype Fitdef ;
  size_t *sims ;
} ;

struct stat_inputs {
  resample_type restype ;
  size_t Nsamples ;
} ;

// input parameters
struct input_params {
  struct fit_info Fit ;
  struct stat_inputs *In_Stats ;
  size_t Ntraj ;
  struct traj *Traj ;
} ;

#endif
