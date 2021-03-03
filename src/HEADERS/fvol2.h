#ifndef FVOL2_H
#define FVOL2_H

int fvol2_NMAX( void ) ;
double
ffvol2( const struct x_desc X , const double *fparams , const size_t Npars ) ;
void
fvol2_f( double *f , const void *data , const double *fparams ) ;
void
fvol2_df( double **df , const void *data , const double *fparams ) ;
void
fvol2_d2f( double **d2f , const void *data , const double *fparams ) ;
void
fvol2_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit ) ;

#endif
