#ifndef SU2_SHITFIT_H
#define SU2_SHITFIT_H

double
fsu2_shitfit( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
su2_shitfit_f( double *f , const void *data , const double *fparams ) ;

// derivatives
void
su2_shitfit_df( double **df , const void *data , const double *fparams ) ;

// second derivatives are all zero
void
su2_shitfit_d2f( double **d2f , const void *data , const double *fparams ) ;

// polynomial guesses
void
su2_shitfit_guesses( double *fparams ,
		     const struct data_info Data ,
		     const struct fit_info Fit ) ;

#endif
