#ifndef FVOL_DELTA_FITT0V2_H
#define FVOL_DELTA_FITT0V2_H

void
init_phi3v2( const size_t Nboots ) ;

void set_phi3v2( const size_t sample_idx , const bool is_avg ) ;

double
ffvol_deltav2( const struct x_desc X , const double *fparams , const size_t Npars ) ;
void
fvol_deltav2_f( double *f , const void *data , const double *fparams ) ;
void
fvol_deltav2_df( double **df , const void *data , const double *fparams ) ;
void
fvol_deltav2_d2f( double **d2f , const void *data , const double *fparams ) ;
void
fvol_deltav2_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit ) ;

#endif
