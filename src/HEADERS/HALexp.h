#ifndef HALEXP_H
#define HALEXP_H

double
fHALexp( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
HALexp_f( double *f , const void *data , const double *fparams ) ;

void
HALexp_df( double **df , const void *data , const double *fparams ) ;

void
HALexp_d2f( double **d2f , const void *data , const double *fparams ) ;

void
HALexp_guesses( double *fparams ,
	     const struct data_info Data ,
	     const struct fit_info Fit ) ;

#endif
