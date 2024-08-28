#ifndef FVOL_DELTA_FITT0_H
#define FVOL_DELTA_FITT0_H

void
init_phi3( const size_t Nboots ) ;
void set_phi3( const size_t sample_idx , const bool is_avg ) ;

double
ffvol_delta( const struct x_desc X , const double *fparams , const size_t Npars ) ;
void
fvol_delta_f( double *f , const void *data , const double *fparams ) ;
void
fvol_delta_df( double **df , const void *data , const double *fparams ) ;
void
fvol_delta_d2f( double **d2f , const void *data , const double *fparams ) ;
void
fvol_delta_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit ) ;

#endif
