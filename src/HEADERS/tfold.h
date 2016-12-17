#ifndef TFOLD_H
#define TFOLD_H

int
time_fold( struct resampled *sample ,
	   const double complex *C ,
	   const size_t LT ,
	   const foldtype fold ,
	   const size_t meas ) ;

#endif
