#ifndef ADLER_ALPHA_D0_H
#define ADLER_ALPHA_D0_H

void
set_mu_adleralpha( const double munew ) ;

void
set_Q1_adleralpha( const double Q1new ) ;

double
fadleralpha_D0( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
adleralpha_D0_f( double *f , const void *data , const double *fparams ) ;

void
adleralpha_D0_df( double **df , const void *data , const double *fparams ) ;

void
adleralpha_D0_d2f( double **d2f , const void *data , const double *fparams ) ;

void
adleralpha_D0_guesses( double *fparams ,
		       const struct data_info Data ,
		       const struct fit_info Fit ) ;
#endif
