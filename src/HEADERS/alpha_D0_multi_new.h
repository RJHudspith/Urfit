#ifndef ALPHA_D0_MULTI_NEW_H
#define ALPHA_D0_MULTI_NEW_H

void
set_Q1_multi2( const double val , const size_t idx ) ;

void
set_mu_multi2( const double munew ) ;

double
falpha_D0_multi2( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
alpha_D0_multi2_f( double *f , const void *data , const double *fparams ) ;

void
alpha_D0_multi2_df( double **df , const void *data , const double *fparams ) ;

void
alpha_D0_multi2_d2f( double **d2f , const void *data , const double *fparams ) ;

void
alpha_D0_multi2_guesses( double *fparams ,
		  const struct data_info Data ,
		  const struct fit_info Fit ) ;
#endif
