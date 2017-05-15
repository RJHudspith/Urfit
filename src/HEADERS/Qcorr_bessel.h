#ifndef QCORR_BESSEL_H
#define QCORR_BESSEL_H

double
fQcorr_bessel( const struct x_desc X , const double *fparams , const size_t Npars ) ;
  
void
Qcorr_bessel_f( double *f , const void *data , const double *fparams ) ;

void
Qcorr_bessel_df( double **df , const void *data , const double *fparams ) ;

void
Qcorr_bessel_d2f( double **d2f , const void *data , const double *fparams ) ;

void
Qcorr_bessel_guesses( double *fparams ,
		      const struct data_info Data ,
		      const struct fit_info Fit ) ;

#endif
