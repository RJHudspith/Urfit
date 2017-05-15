#ifndef HISTOGRAM_H
#define HISTOGRAM_H

struct histogram*
histogram( const double *data ,
	   const size_t Ndata ,
	   const size_t Nbins ) ;

double
hist_min( size_t *min ,
	  const struct histogram *hist ,
	  const size_t Nbin ) ;

#endif
