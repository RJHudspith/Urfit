#ifndef FITSOL2_H
#define FITSOL2_H

void
set_psq_sol2( const double *p2 ,
	     const size_t NP2 ) ;

double
fsol2( const struct x_desc X ,
      const double *fparams ,
      const size_t Npars ) ;

void
sol2_f( double *f ,
	      const void *data ,
	      const double *fparams ) ;

void
sol2_df( double **df ,
	       const void *data ,
	       const double *fparams ) ;
void
sol2_d2f( double **d2f ,
		const void *data ,
		const double *fparams ) ;

void
sol2_guesses( double *fparams ,
		    const struct data_info Data ,
		    const struct fit_info Fit ) ;

#endif
