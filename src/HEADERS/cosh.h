#ifndef COSH_H
#define COSH_H

double
fcosh( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
cosh_f( double *f , const void *data , const double *fparams ) ;

void
cosh_df( double **df , const void *data , const double *fparams ) ;

void
cosh_d2f( double **d2f , const void *data , const double *fparams ) ;

void
cosh_guesses( double *fparams ) ;

void
cosh_priors( double *priors , double *err_priors ) ;

#endif
