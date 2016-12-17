#ifndef POLY_H
#define POLY_H

double
fpoly( const struct x_desc X , const double *fparams ) ;

void
poly_f( double *f , const void *data , const double *fparams ) ;

void
poly_df( double **df , const void *data , const double *fparams ) ;

void
poly_d2f( double **d2f , const void *data , const double *fparams ) ;

void
poly_guesses( double *fparams ) ;

void
poly_priors( double *priors , double *err_priors ) ;

#endif
