#ifndef ZV_EXP_H
#define ZV_EXP_H

double
fZV_exp( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
ZV_exp_f( double *f , const void *data , const double *fparams ) ;

void
ZV_exp_df( double **df , const void *data , const double *fparams ) ;

void
ZV_exp_d2f( double **d2f , const void *data , const double *fparams ) ;

void
ZV_exp_guesses( double *fparams ,
		const struct data_info Data ,
		const struct fit_info Fit ) ;

#endif
