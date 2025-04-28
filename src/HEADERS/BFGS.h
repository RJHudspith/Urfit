#ifndef BFGS_H
#define BFGS_H

int
BFGS_iter( void *fdesc ,
	   const void *data ,
	   const double **W ,
	   const double TOL ) ;

#endif
