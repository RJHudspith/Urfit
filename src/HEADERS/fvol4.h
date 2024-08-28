#ifndef FVOL4_H
#define FVOL4_H

double
ffvol4( const struct x_desc X , const double *fparams , const size_t Npars ) ;
void
fvol4_f( double *f , const void *data , const double *fparams ) ;
void
fvol4_df( double **df , const void *data , const double *fparams ) ;
void
fvol4_d2f( double **d2f , const void *data , const double *fparams ) ;
void
fvol4_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit ) ;

#endif
