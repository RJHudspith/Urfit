#ifndef LINE_SEARCH_H
#define LINE_SEARCH_H

// NRC min function for Powell's
double
linmin( const int n ,
	double p[n] ,
	double xi[n] ,
	double *fret ,
	struct ffunction *f2 , 
	const struct fit_descriptor *fdesc ,
	const double **W ,
	const void *data ) ;

// my rougher backtrack+golden code
double
line_search( struct ffunction *f2 ,
	     const struct ffunction f1 ,
	     const double *grad ,
	     const double *descent ,
	     const struct fit_descriptor fdesc ,
	     const void *data ,
	     const double **W ) ;

// gets the derivative of the \chi^2 function
void
get_gradient( double *grad ,
	      const double **W ,
	      const struct fit_descriptor *Fit ) ;
  
#endif
