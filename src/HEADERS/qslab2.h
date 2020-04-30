#ifndef QSLAB2_H
#define QSLAB2_H

double
fqslab2( const struct x_desc X , const double *fparams , const size_t Npars ) ;
  
void
qslab2_f( double *f , const void *data , const double *fparams ) ;

void
qslab2_df( double **df , const void *data , const double *fparams ) ;

void
qslab2_d2f( double **d2f , const void *data , const double *fparams ) ;

void
qslab2_guesses( double *fparams ,
	     const struct data_info Data ,
	       const struct fit_info Fit ) ;
#endif
