#ifndef SUN_CONT_H
#define SUN_CONT_H
double
fSUN_cont( const struct x_desc X , const double *fparams , const size_t Npars );
void
SUN_cont_f( double *f , const void *data , const double *fparams );
void
SUN_cont_df( double **df , const void *data , const double *fparams );
void
SUN_cont_d2f( double **d2f , const void *data , const double *fparams );
void
SUN_cont_guesses( double *fparams ,
		  const struct data_info Data ,
		  const struct fit_info Fit ) ;
#endif
