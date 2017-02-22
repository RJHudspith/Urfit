#ifndef POLY_H
#define POLY_H

double
fpoly( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
poly_f( double *f , const void *data , const double *fparams ) ;

void
poly_df( double **df , const void *data , const double *fparams ) ;

void
poly_d2f( double **d2f , const void *data , const double *fparams ) ;

void
poly_guesses( double *fparams ,
	      const struct data_info Data ,
	      const struct fit_info Fit ) ;

void
poly_linmat( double **U ,
	     const void *data ,
	     const size_t Nparam ,
	     const size_t Nlogic ) ;

#endif
