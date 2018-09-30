#ifndef C4C7_H
#define C4C7_H
double
fc4c7( const struct x_desc X , const double *fparams , const size_t Npars );
void
c4c7_f( double *f , const void *data , const double *fparams );
void
c4c7_df( double **df , const void *data , const double *fparams );
void
c4c7_d2f( double **d2f , const void *data , const double *fparams );
void
c4c7_guesses( double *fparams ,
	      const struct data_info Data ,
	      const struct fit_info Fit ) ;
#endif
