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

struct resampled *
resample_data( const struct resampled *RAW ,
	       const size_t N ,
	       const resample_type restype ,
	       const int NBOOTS ) ;

struct resampled **
resample_array( const struct resampled **RAW ,
		const size_t NSLICES ,
		const size_t *NDATA ,
		const size_t resample , 
		const size_t NBOOTS ) ;

#endif
