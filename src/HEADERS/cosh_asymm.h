#ifndef COSH_ASYMM_H
#define COSH_ASYMM_H

double
fcosh_asymm( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
cosh_asymm_f( double *f , const void *data , const double *fparams ) ;

void
cosh_asymm_df( double **df , const void *data , const double *fparams ) ;

void
cosh_asymm_d2f( double **d2f , const void *data , const double *fparams ) ;

#endif
