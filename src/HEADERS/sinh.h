#ifndef SINH_H
#define SINH_H

double
fsinh( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
sinh_f( double *f , const void *data , const double *fparams ) ;

void
sinh_df( double **df , const void *data , const double *fparams ) ;

void
sinh_d2f( double **d2f , const void *data , const double *fparams ) ;

void
sinh_guesses( double *fparams ) ;

void
sinh_priors( double *priors , double *err_priors ) ;

#endif
