#ifndef COSH_PLUSC_H
#define COSH_PLUSC_H

double
fcosh_plusc( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
cosh_plusc_f( double *f , const void *data , const double *fparams ) ;

void
cosh_plusc_df( double **df , const void *data , const double *fparams ) ;

void
cosh_plusc_d2f( double **d2f , const void *data , const double *fparams ) ;

void
cosh_plusc_guesses( double *fparams , const size_t Nlogic ) ;

#endif
