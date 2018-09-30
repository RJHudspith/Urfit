#ifndef ADLER_ALPHA_D0_MULTI_H
#define ADLER_ALPHA_D0_MULTI_H

void
set_mu_multi_adler( const double munew ) ;

double
fadleralpha_D0_multi( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
adleralpha_D0_multi_f( double *f , const void *data , const double *fparams ) ;

void
adleralpha_D0_multi_df( double **df , const void *data , const double *fparams ) ;

void
adleralpha_D0_multi_d2f( double **d2f , const void *data , const double *fparams ) ;

void
adleralpha_D0_multi_guesses( double *fparams ,
			     const struct data_info Data ,
			     const struct fit_info Fit ) ;
#endif
