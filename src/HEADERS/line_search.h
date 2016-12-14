#ifndef LINE_SEARCH_H
#define LINE_SEARCH_H

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

double
line_search2( struct ffunction *f2 ,
	     const struct ffunction f1 ,
	     const double *grad ,
	     const double *descent ,
	     const struct fit_descriptor fdesc ,
	     const void *data ,
	     const double **W ,
	     const double alpha ) ;

#endif
