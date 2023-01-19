#ifndef CORNELLV2_H
#define CORNELLV2_H

double
fcornellv2( const struct x_desc X , const double *fparams , const size_t Npars ) ;
void
cornellv2_f( double *f , const void *data , const double *fparams ) ;
void
cornellv2_df( double **df , const void *data , const double *fparams ) ;
void
cornellv2_d2f( double **d2f , const void *data , const double *fparams ) ;
void
cornellv2_guesses( double *fparams ,
	       const struct data_info Data ,
		 const struct fit_info Fit ) ;
#endif
