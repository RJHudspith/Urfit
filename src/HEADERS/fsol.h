#ifndef FITSOL_H
#define FITSOL_H

void
set_psq_sol( const double *p2 ,
	     const size_t NP2 ) ;

double
fsol( const struct x_desc X ,
      const double *fparams ,
      const size_t Npars ) ;

void
sol_f( double *f ,
	      const void *data ,
	      const double *fparams ) ;

void
sol_df( double **df ,
	       const void *data ,
	       const double *fparams ) ;
void
sol_d2f( double **d2f ,
		const void *data ,
		const double *fparams ) ;

void
sol_guesses( double *fparams ,
		    const struct data_info Data ,
		    const struct fit_info Fit ) ;

#endif
