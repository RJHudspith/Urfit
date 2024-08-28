#ifndef FVOL5_H
#define FVOL5_H

double
ffvol5( const struct x_desc X , const double *fparams , const size_t Npars ) ;
void
fvol5_f( double *f , const void *data , const double *fparams ) ;
void
fvol5_df( double **df , const void *data , const double *fparams ) ;
void
fvol5_d2f( double **d2f , const void *data , const double *fparams ) ;
void
fvol5_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit ) ;

#endif
