#ifndef GENS_H
#define GENS_H

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>

typedef enum { 
  UNWEIGHTED , UNCORRELATED , CORRELATED
} corrtype ;

struct data {
  size_t n;
  double *y;
};

struct ffunction {
  double *f ; // ( y_i - f_i(a) )
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
  void (*F) ( double *f , const void *data , const double *fparams ) ;
  void (*dF) ( double **df , const void *data , const double *fparams ) ;
  void (*d2F) ( double **d2f , const void *data , const double *fparams ) ;
  void (*guesses) ( double *fparams ) ;
  void (*set_priors) ( double *priors , double *err_priors ) ;
  size_t NPARAMS ;
} ;

#define UNINIT_FLAG (123456789)

#endif
