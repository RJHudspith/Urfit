#ifndef TEST_H
#define TEST_H

double
ftest( const struct x_desc X , const double *fparams , const size_t Npars ) ;

void
test_f( double *f , const void *data , const double *fparams ) ;

void
test_df( double **df , const void *data , const double *fparams ) ;

void
test_d2f( double **d2f , const void *data , const double *fparams ) ;

void
test_guesses( double *fparams ,
	       const struct data_info Data ,
	      const struct fit_info Fit ) ;

#endif
