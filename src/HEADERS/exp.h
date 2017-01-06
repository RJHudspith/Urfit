#ifndef EXP_H
#define EXP_H

double
fexp( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
exp_f( double *f , const void *data , const double *fparams ) ;

void
exp_df( double **df , const void *data , const double *fparams ) ;

void
exp_d2f( double **d2f , const void *data , const double *fparams ) ;

void
exp_guesses( double *fparams , const size_t Nlogic ) ;

void
exp_priors( double *priors , double *err_priors ) ;

#endif
