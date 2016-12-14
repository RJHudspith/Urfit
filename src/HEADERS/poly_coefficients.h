#ifndef POLY_COEFFICIENTS_H
#define POLY_COEFFICIENTS_H

int
compute_coefficients( double *coeffs ,
		      double *chisq ,
		      const double *y ,
		      const double *sigma ,
		      const double *x ,
		      const size_t M ,
		      const size_t N ) ;

void
write_polynomial( const double *__restrict coeffs ,
		  const size_t POLY_ORDER ) ;

#endif 
