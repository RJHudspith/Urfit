#ifndef POWELL_H
#define POWELL_H

int
powell_iter( void *fdesc ,
	     const void *data ,
	     const double **W ,
	     const double TOL ) ;

#endif
