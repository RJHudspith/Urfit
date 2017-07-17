#ifndef PPAA_H
#define PPAA_H


double
fppaa( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
ppaa_f( double *f , const void *data , const double *fparams ) ;

// derivatives
void
ppaa_df( double **df , const void *data , const double *fparams ) ;

void
ppaa_d2f( double **d2f , const void *data , const double *fparams ) ;

void
ppaa_guesses( double *fparams ,
	      const struct data_info Data ,
	      const struct fit_info Fit ) ;

#endif
