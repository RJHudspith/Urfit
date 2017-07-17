#ifndef QSLAB_H
#define QSLAB_H

double
fqslab( const struct x_desc X , const double *fparams , const size_t Npars ) ;
  
void
qslab_f( double *f , const void *data , const double *fparams ) ;

void
qslab_df( double **df , const void *data , const double *fparams ) ;

void
qslab_d2f( double **d2f , const void *data , const double *fparams ) ;

void
qslab_guesses( double *fparams ,
	     const struct data_info Data ,
	       const struct fit_info Fit ) ;
#endif
