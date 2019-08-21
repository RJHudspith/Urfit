#ifndef PEXP_H
#define PEXP_H

double
fPexp( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
Pexp_f( double *f , const void *data , const double *fparams ) ;

// derivatives
void
Pexp_df( double **df , const void *data , const double *fparams ) ;

// second derivatives
void
Pexp_d2f( double **d2f , const void *data , const double *fparams ) ;

void
Pexp_guesses( double *fparams , const size_t Nlogic ) ;

#endif
