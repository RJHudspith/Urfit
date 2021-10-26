#ifndef FVOLCC_H
#define FVOLCC_H

double
ffvolcc( const struct x_desc X , const double *fparams , const size_t Npars ) ;
void
fvolcc_f( double *f , const void *data , const double *fparams ) ;
void
fvolcc_df( double **df , const void *data , const double *fparams ) ;
void
fvolcc_d2f( double **d2f , const void *data , const double *fparams ) ;
void
fvolcc_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit ) ;

#endif
