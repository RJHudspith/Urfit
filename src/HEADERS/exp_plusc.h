#ifndef EXP_PLUSC_H
#define EXP_PLUSC_H

double
fexp_plusc( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
exp_plusc_f( double *f , const void *data , const double *fparams ) ;

void
exp_plusc_df( double **df , const void *data , const double *fparams ) ;

void
exp_plusc_d2f( double **d2f , const void *data , const double *fparams ) ;

void
exp_plusc_guesses( double *fparams , const size_t Nlogic ) ;

void
exp_plusc_priors( double *priors , double *err_priors ) ;

#endif
