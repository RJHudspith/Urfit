#ifndef ALPHA_D0_multi_H
#define ALPHA_D0_multi_H

void
set_Q1_multi( const double val , const size_t idx ) ;

void
set_mu_multi( const double munew ) ;

double
falpha_D0_multi( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
alpha_D0_multi_f( double *f , const void *data , const double *fparams ) ;

void
alpha_D0_multi_df( double **df , const void *data , const double *fparams ) ;

void
alpha_D0_multi_d2f( double **d2f , const void *data , const double *fparams ) ;

void
alpha_D0_multi_guesses( double *fparams ,
		  const struct data_info Data ,
		  const struct fit_info Fit ) ;
#endif
