#ifndef RAW_H
#define RAW_H

void
average( double *ave , double *err , 
	 const double *meas , const size_t N ) ;

void
raw_err( struct resampled *replicas ) ;

#endif
