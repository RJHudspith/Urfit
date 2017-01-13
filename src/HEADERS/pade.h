#ifndef PADE_H
#define PADE_H

void 
pade_set_nm( const size_t new_n ,
	     const size_t new_m ) ;

void
pade_get_nm( size_t *new_n ,
	     size_t *new_m ) ;

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
pade_guesses( double *fparams , const size_t Nlogic ) ;

#endif
