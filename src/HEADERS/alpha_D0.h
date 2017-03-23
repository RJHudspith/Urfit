#ifndef ALPHA_D0_H
#define ALPHA_D0_H

void
set_Q1( const double val , const size_t idx ) ;

double
falpha_D0( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
alpha_D0_f( double *f , const void *data , const double *fparams ) ;

void
alpha_D0_df( double **df , const void *data , const double *fparams ) ;

void
alpha_D0_d2f( double **d2f , const void *data , const double *fparams ) ;

void
alpha_D0_guesses( double *fparams ,
		  const struct data_info Data ,
		  const struct fit_info Fit ) ;
#endif
