#ifndef PADE_H
#define PADE_H

double
fpade( const struct x_desc X , 
       const double *fparams ,
       const size_t Npars ) ;

void
pade_f( double *f , 
	const void *data , 
	const double *fparams ) ;

void
pade_df( double **df , 
	 const void *data , 
	 const double *fparams ) ;

void
pade_d2f( double **d2f , 
	  const void *data , 
	  const double *fparams ) ;

void
pade_guesses( double *fparams ,
	      const struct data_info Data ,
	      const struct fit_info Fit ) ;

#endif
