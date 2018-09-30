#ifndef FVOL1_H
#define FVOL1_H

double
ffvol1( const struct x_desc X , const double *fparams , const size_t Npars ) ;
void
fvol1_f( double *f , const void *data , const double *fparams ) ;
void
fvol1_df( double **df , const void *data , const double *fparams ) ;
void
fvol1_d2f( double **d2f , const void *data , const double *fparams ) ;
void
fvol1_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit ) ;

#endif
