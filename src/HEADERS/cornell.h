#ifndef CORNELL_H
#define CORNELL_H

double
fcornell( const struct x_desc X , const double *fparams , const size_t Npars ) ;
void
cornell_f( double *f , const void *data , const double *fparams ) ;
void
cornell_df( double *f , const void *data , const double *fparams ) ;
void
cornell_d2f( double **d2f , const void *data , const double *fparams ) ;
void
cornell_guesses( double *fparams ,
	       const struct data_info Data ,
		 const struct fit_info Fit ) ;
#endif
