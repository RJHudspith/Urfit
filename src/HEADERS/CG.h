#ifndef CG_H
#define CG_H

int
cg_iter( struct fit_descriptor *fdesc ,
	 const void *data ,
	 const double **W ,
	 const double TOL ) ;

#endif
