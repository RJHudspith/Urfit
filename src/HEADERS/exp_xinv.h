#ifndef EXP_XINV_H
#define EXP_XINV_H

double
fexp_xinv( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
exp_xinv_f( double *f , const void *data , const double *fparams ) ;

void
exp_xinv_df( double **df , const void *data , const double *fparams ) ;

void
exp_xinv_d2f( double **d2f , const void *data , const double *fparams ) ;

void
exp_xinv_guesses( double *fparams ,
	     const struct data_info Data ,
	     const struct fit_info Fit ) ;

#endif
