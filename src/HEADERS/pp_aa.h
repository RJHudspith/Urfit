#ifndef PP_AA_H
#define PP_AA_H


double
fpp_aa( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
pp_aa_f( double *f , const void *data , const double *fparams ) ;

// derivatives
void
pp_aa_df( double **df , const void *data , const double *fparams ) ;

void
pp_aa_d2f( double **d2f , const void *data , const double *fparams ) ;

void
pp_aa_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit ) ;

#endif
