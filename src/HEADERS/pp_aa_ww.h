#ifndef PP_AA_WW_H
#define PP_AA_WW_H

double
fpp_aa_ww( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
pp_aa_ww_f( double *f , const void *data , const double *fparams ) ;

// derivatives
void
pp_aa_ww_df( double **df , const void *data , const double *fparams ) ;

void
pp_aa_ww_d2f( double **d2f , const void *data , const double *fparams ) ;

void
pp_aa_ww_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit ) ;

#endif
