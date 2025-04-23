#ifndef LINE_SEARCH_H
#define LINE_SEARCH_H

void
linmin( const int n ,
	double p[n] ,
	double xi[n] ,
	double *fret ,
	struct ffunction *f2 , 
	const struct fit_descriptor *fdesc ,
	const double **W ,
	const void *data ) ;

double
line_search( struct ffunction *f2 ,
	     const struct ffunction f1 ,
	     const double *grad ,
	     const double *descent ,
	     const struct fit_descriptor fdesc ,
	     const void *data ,
	     const double **W ,
	     const size_t jidx ,
	     const double alpha ) ;

#endif
