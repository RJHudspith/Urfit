#ifndef POLES_H
#define POLES_H

double
fpoles( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
poles_f( double *f , const void *data , const double *fparams ) ;

// derivatives
void
poles_df( double **df , const void *data , const double *fparams ) ;

void
poles_d2f( double **d2f , const void *data , const double *fparams ) ;

void
poles_linmat( double **U ,
	      const void *data ,
	      const size_t N ,
	      const size_t M ,
	      const size_t Nlogic ) ;

void
poles_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit ) ;
#endif
