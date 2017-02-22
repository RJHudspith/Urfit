#ifndef GLS_H
#define GLS_H

int
gls_iter( void *fdesc ,
	  const void *data ,
	  const double **W ,
	  const double TOL ) ;

#endif
