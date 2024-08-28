#ifndef FVOL6_H
#define FVOL6_H

double
ffvol6( const struct x_desc X , const double *fparams , const size_t Npars ) ;
void
fvol6_f( double *f , const void *data , const double *fparams ) ;
void
fvol6_df( double **df , const void *data , const double *fparams ) ;
void
fvol6_d2f( double **d2f , const void *data , const double *fparams ) ;
void
fvol6_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit ) ;

#endif
