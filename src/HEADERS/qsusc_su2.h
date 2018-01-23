#ifndef QSUSC_SU2_H
#define QSUSC_SU2_H

double
fqsusc_su2( const struct x_desc X , const double *fparams , const size_t Npars ) ;
void
qsusc_su2_f( double *f , const void *data , const double *fparams ) ;
void
qsusc_su2_df( double **df , const void *data , const double *fparams ) ;
void
qsusc_su2_d2f( double **d2f , const void *data , const double *fparams ) ;
void
qsusc_su2_guesses( double *fparams ,
		   const struct data_info Data ,
		   const struct fit_info Fit ) ;

#endif
