#ifndef POLY_H
#define POLY_H

void
poly_set_n( const size_t new_n ) ;

double
fpoly( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
poly_f( double *f , const void *data , const double *fparams ) ;

void
poly_df( double **df , const void *data , const double *fparams ) ;

void
poly_d2f( double **d2f , const void *data , const double *fparams ) ;

void
poly_guesses( double *fparams , const size_t Nlogic ) ;

#endif
