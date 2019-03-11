#ifndef TANH_H
#define TANH_H

double
ftanh( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
tanh_f( double *f , const void *data , const double *fparams ) ;

void
tanh_df( double **df , const void *data , const double *fparams ) ;

void
tanh_d2f( double **d2f , const void *data , const double *fparams ) ;

#endif
