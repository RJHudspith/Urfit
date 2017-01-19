#ifndef SUMMATION_H
#define SUMMATION_H

double
kahan_summation( const double *data ,
		 const size_t Ndata ) ;

double
knuth_average( double *err ,
	       const double *data ,
	       const size_t Ndata ) ;

#endif
