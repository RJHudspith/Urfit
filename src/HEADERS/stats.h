#ifndef STATS_H
#define STATS_H

void
average( double *ave , 
	 double *err , 
	 const double *meas , 
	 const size_t N ) ;

struct resampled
bin_the_data( const struct resampled RAW ,
	      const size_t binning ) ;

void
compute_err( struct resampled *replicas ) ;

int
resample_data( struct input_params *Input ) ;

#endif
