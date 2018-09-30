#ifndef UDCB_HEAVY_H
#define UDCB_HEAVY_H

double
fudcb_heavy( const struct x_desc X , const double *fparams , const size_t Npars ) ;
void
udcb_heavy_f( double *f , const void *data , const double *fparams ) ;
void
udcb_heavy_df( double **df , const void *data , const double *fparams ) ;
void
udcb_heavy_d2f( double **d2f , const void *data , const double *fparams ) ;
void
udcb_heavy_guesses( double *fparams ,
		    const struct data_info Data ,
		    const struct fit_info Fit ) ;

#endif
