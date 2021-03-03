#ifndef FVOL3_H
#define FVOL3_H

double
ffvol3( const struct x_desc X , const double *fparams , const size_t Npars ) ;
void
fvol3_f( double *f , const void *data , const double *fparams ) ;
void
fvol3_df( double **df , const void *data , const double *fparams ) ;
void
fvol3_d2f( double **d2f , const void *data , const double *fparams ) ;
void
fvol3_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit ) ;

#endif
