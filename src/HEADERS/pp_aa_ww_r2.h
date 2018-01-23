#ifndef PP_AA_WW_R2_H
#define PP_AA_WW_R2_H

double
fpp_aa_ww_r2( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
pp_aa_ww_r2_f( double *f , const void *data , const double *fparams ) ;

// derivatives
void
pp_aa_ww_r2_df( double **df , const void *data , const double *fparams ) ;

void
pp_aa_ww_r2_d2f( double **d2f , const void *data , const double *fparams ) ;

void
pp_aa_ww_r2_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit ) ;

#endif
