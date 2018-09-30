#ifndef NRQCD_EXP_H
#define NRQCD_EXP_H

void
set_psq( const double *p2 ,
	 const size_t NP2 ) ;

double
fnrqcd_exp( const struct x_desc X ,
	    const double *fparams ,
	    const size_t Npars ) ;

void
nrqcd_exp_f( double *f ,
	     const void *data ,
	     const double *fparams ) ;

void
nrqcd_exp_df( double **df ,
	      const void *data ,
	      const double *fparams ) ;
void
nrqcd_exp_d2f( double **d2f , const void *data , const double *fparams ) ;

void
nrqcd_exp_guesses( double *fparams ,
		   const struct data_info Data ,
		   const struct fit_info Fit ) ;

#endif
