#ifndef PP_AA_EXP_H
#define PP_AA_EXP_H

double
fpp_aa_exp( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
pp_aa_exp_f( double *f , const void *data , const double *fparams ) ;

// derivatives
void
pp_aa_exp_df( double **df , const void *data , const double *fparams ) ;

void
pp_aa_exp_d2f( double **d2f , const void *data , const double *fparams ) ;

void
pp_aa_exp_guesses( double *fparams ,
	       const struct data_info Data ,
	       const struct fit_info Fit ) ;

#endif
