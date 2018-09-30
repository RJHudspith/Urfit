#ifndef NRQCD_EXP2_H
#define NRQCD_EXP2_H

void
set_psq_nrqcd2( const double *p2 ,
		const size_t NP2 ) ;

double
fnrqcd_exp2( const struct x_desc X ,
	     const double *fparams ,
	     const size_t Npars ) ;

void
nrqcd_exp2_f( double *f ,
	      const void *data ,
	      const double *fparams ) ;

void
nrqcd_exp2_df( double **df ,
	       const void *data ,
	       const double *fparams ) ;
void
nrqcd_exp2_d2f( double **d2f ,
		const void *data ,
		const double *fparams ) ;

void
nrqcd_exp2_guesses( double *fparams ,
		    const struct data_info Data ,
		    const struct fit_info Fit ) ;

#endif
